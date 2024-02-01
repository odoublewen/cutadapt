import errno
import io
import sys
from abc import ABC
from enum import Enum
from typing import BinaryIO, Optional, Dict, Tuple, List, TextIO, Union

import dnaio
from xopen import xopen

from cutadapt.utils import logger

try:
    import resource
except ImportError:
    # Windows
    resource = None  # type: ignore


def xopen_rb_raise_limit(path: str):
    """
    Open a (possibly compressed) file for reading in binary mode, trying to avoid the
    "Too many open files" problem using `open_raise_limit`.
    """
    mode = "rb"
    f = open_raise_limit(xopen, path, mode, threads=0)
    logger.debug("Opening '%s', mode '%s' with xopen resulted in %s", path, mode, f)
    return f


def open_raise_limit(func, *args, **kwargs):
    """
    Run 'func' (which should be some kind of open() function) and return its result.
    If "Too many open files" occurs, increase limit and try again.
    """
    try:
        f = func(*args, **kwargs)
    except OSError as e:
        if e.errno == errno.EMFILE:  # Too many open files
            logger.debug("Too many open files, attempting to raise soft limit")
            raise_open_files_limit(8)
            f = func(*args, **kwargs)
        else:
            raise
    return f


def raise_open_files_limit(n):
    if resource is None:
        return
    soft, hard = resource.getrlimit(resource.RLIMIT_NOFILE)
    soft = min(soft + n, hard)
    resource.setrlimit(resource.RLIMIT_NOFILE, (soft, hard))


class FileOpener:
    def __init__(self, compression_level: int = 6, threads: Optional[int] = None):
        """
        threads -- no. of external compression threads.
            0: write in-process
            None: min(cpu_count(), 4)
        """
        self.compression_level = compression_level
        self.threads = threads

    def xopen(self, path, mode):
        threads = self.threads if "w" in mode else 0
        f = open_raise_limit(
            xopen, path, mode, compresslevel=self.compression_level, threads=threads
        )
        if "w" in mode:
            extra = f" (compression level {self.compression_level}, {threads} threads)"
        else:
            extra = ""
        logger.debug(
            "Opening '%s', mode '%s'%s with xopen resulted in %s", path, mode, extra, f
        )
        return f

    def xopen_or_none(self, path, mode):
        """Return opened file or None if the path is None"""
        if path is None:
            return None
        return self.xopen(path, mode)

    def xopen_pair(self, path1: str, path2: Optional[str], mode):
        if path1 is None and path2 is not None:
            raise ValueError(
                "When giving paths for paired-end files, only providing the second"
                " file is not supported"
            )
        file1 = self.xopen_or_none(path1, mode)
        file2 = self.xopen_or_none(path2, mode)
        return file1, file2

    def dnaio_open(self, *args, **kwargs):
        kwargs["opener"] = self.xopen
        f = dnaio.open(*args, **kwargs)
        if not isinstance(args[0], io.BytesIO):
            logger.debug(
                "Opening %r, mode '%s' with dnaio resulted in %s",
                args[0],
                kwargs["mode"],
                f,
            )
        return f

    def dnaio_open_raise_limit(self, *args, **kwargs):
        """
        Open a FASTA/FASTQ file for writing. If it fails because the number of open files
        would be exceeded, try to raise the soft limit and re-try.
        """
        return open_raise_limit(self.dnaio_open, *args, **kwargs)


class InputFiles:
    def __init__(
        self,
        *files: BinaryIO,
        interleaved: bool = False,
    ):
        self._files = files
        self.interleaved = interleaved
        for f in self._files:
            assert f is not None

    def open(self):
        return dnaio.open(*self._files, interleaved=self.interleaved, mode="r")

    def close(self) -> None:
        for file in self._files:
            file.close()


class InputPaths:
    def __init__(self, *paths: str, interleaved: bool = False):
        self.paths = paths
        self.interleaved = interleaved

    def open(self) -> InputFiles:
        files = [xopen_rb_raise_limit(path) for path in self.paths]
        return InputFiles(*files, interleaved=self.interleaved)


class ProxyWriter(ABC):
    def drain(self) -> List[bytes]:
        pass


class ProxyTextFile(ProxyWriter):
    def __init__(self):
        self._buffer = io.BytesIO()
        self._file = io.TextIOWrapper(self._buffer)

    def write(self, text):
        self._file.write(text)

    def drain(self) -> List[bytes]:
        self._file.flush()
        chunk = self._buffer.getvalue()
        self._buffer.seek(0)
        self._buffer.truncate()
        return [chunk]

    def __getstate__(self):
        """TextIOWrapper cannot be pickled. Just don’t include our state."""
        return True  # ensure __setstate__ is called

    def __setstate__(self, state):
        self.__init__()


class ProxyRecordWriter(ProxyWriter):
    def __init__(self, n_files: int, **kwargs):
        self._n_files = n_files
        self._kwargs = kwargs
        self._buffers = [io.BytesIO() for _ in range(n_files)]
        self._writer = open_raise_limit(dnaio.open, *self._buffers, mode="w", **kwargs)

    def write(self, *args, **kwargs):
        self._writer.write(*args, **kwargs)

    def drain(self) -> List[bytes]:
        chunks = [buf.getvalue() for buf in self._buffers]
        for buf in self._buffers:
            buf.seek(0)
            buf.truncate()
        return chunks

    def __getstate__(self):
        """Exclude the dnaio Reader class from the state"""
        return (self._n_files, self._kwargs)

    def __setstate__(self, state):
        n_files, kwargs = state
        self.__init__(n_files, **kwargs)


class OutputFiles:
    def __init__(
        self,
        *,
        file_opener: FileOpener,
        proxied: bool,
        qualities: bool,
        interleaved: bool,
    ):
        self._file_opener = file_opener
        # TODO do these actually have to be dicts?
        self._binary_files: Dict[str, BinaryIO] = {}
        self._text_files: Dict[str, TextIO] = {}
        self._writers: Dict = {}
        self._proxy_files: List[ProxyWriter] = []
        self._proxied = proxied
        self._to_close = []
        self._qualities = qualities
        self._interleaved = interleaved

    def open_text(self, path):
        assert path not in self._binary_files  # TODO
        # TODO
        # - serial runner needs only text_file
        # - parallel runner needs binary_file and proxy_file
        # split into SerialOutputFiles and ParallelOutputFiles?
        if self._proxied:
            binary_file = self._file_opener.xopen(path, "wb")
            self._binary_files[path] = binary_file
            proxy_file = ProxyTextFile()
            self._proxy_files.append(proxy_file)
            return proxy_file
        else:
            text_file = self._file_opener.xopen(path, "wt")
            self._text_files[path] = text_file
            return text_file

    def open_record_writer(self, *paths):
        kwargs = dict(qualities=self._qualities)
        if len(paths) not in (1, 2):
            raise ValueError("Expected one or two paths")
        if len(paths) == 2 and paths[1] is None:
            paths = paths[:1]
            kwargs["interleaved"] = True
        for path in paths:
            assert path is not None
            assert path not in self._binary_files  # TODO
        binary_files = []
        for path in paths:
            binary_file = self._file_opener.xopen(path, "wb")
            binary_files.append(binary_file)
            self._binary_files[path] = binary_file
        if self._proxied:
            proxy_writer = ProxyRecordWriter(len(paths), **kwargs)
            self._proxy_files.append(proxy_writer)
            return proxy_writer
        else:
            writer = self._file_opener.dnaio_open(*binary_files, mode="w", **kwargs)
            self._writers[paths] = writer
            return writer

    def open_record_writer_from_binary_io(self, file: BinaryIO, interleaved: bool):
        self._binary_files["fake\0path"] = file  # TODO
        if self._proxied:
            proxy_writer = ProxyRecordWriter(1, qualities=self._qualities, interleaved=interleaved)
            self._proxy_files.append(proxy_writer)
            return proxy_writer
        else:
            writer = self._file_opener.dnaio_open(file, mode="w", qualities=self._qualities, interleaved=interleaved)
            self._writers["fake\0path"] = writer
            return writer

    def binary_files(self):
        return list(self._binary_files.values())

    def proxy_files(self) -> List[ProxyWriter]:
        return self._proxy_files

    def __iter__(self):
        yield from self._binary_files.values()

    def close(self) -> None:
        """Close all output files that are not stdout"""
        # TODO ... that *are not* stdout
        if not self._proxied:
            for f in self._text_files.values():
                f.close()
            for f in self._writers.values():
                f.close()
        for f in self._binary_files.values():
            if f is not sys.stdout.buffer:
                f.close()


class FileFormat(Enum):
    FASTA = 1
    FASTQ = 2
    BAM = 3

    def has_qualities(self) -> bool:
        return self is FileFormat.FASTQ or self is FileFormat.BAM  # TODO BAM?


# TODO copied and adjusted from dnaio; upstream this
def detect_file_format(file: BinaryIO) -> Optional[FileFormat]:
    if file.seekable():
        original_position = file.tell()
        magic = file.read(4)
        file.seek(original_position)
    else:
        # We cannot always use peek() because BytesIO objects do not suppert it
        magic = file.peek(4)[0:4]  # type: ignore
    if magic.startswith(b"@") or magic == b"":
        # Pretend FASTQ for empty input
        return FileFormat.FASTQ
    elif magic.startswith(b">") or magic.startswith(b"#"):
        # Some FASTA variants allow comments
        return FileFormat.FASTA
    elif magic == b"BAM\1":
        return FileFormat.BAM
    return None

