from .histrader import main
from importlib.metadata import version, PackageNotFoundError
try:
    __version__ = version("histrader")
except PackageNotFoundError:
    __version__ = "unknown"