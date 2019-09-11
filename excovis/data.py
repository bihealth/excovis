"""Code for the definition of the data and meta data and loading them from disk.

The types of this modules are used commonly in the excovis app but the actual ``load_`` functions
should only be accessed through the ``.store`` module such that they can be cached by the app.
This module is unaware of the Dash app.
"""

import hashlib
from urllib.parse import urlunparse as _urlunparse
import re

import attr
import fs.path
import fs.tools
from fs.osfs import OSFS
import pysam

from . import settings
from .exceptions import ExcovisException

#: Identifier for fake data
FAKE_DATA_ID = "builtin-fake-data"


@attr.s(auto_attribs=True)
class MetaData:
    """Class to bundle the metadata from ."""

    #: ID (= file name) of the dataset
    id: str
    #: Full path to the file
    path: str
    #: Name of the sample in the BAM file
    sample: str


def redacted_urlunparse(url, redact_with="***"):
    """``urlunparse()`` but redact password."""
    netloc = []
    if url.username:
        netloc.append(url.username)
    if url.password:
        netloc.append(":")
        netloc.append(redact_with)
    if url.hostname:
        if netloc:
            netloc.append("@")
        netloc.append(url.hostname)
    url = url._replace(netloc="".join(netloc))
    return _urlunparse(url)


#: Schemes supported through PyFilesystem
PYFS_SCHEMES = ("file",)


def make_osfs(url):
    """Construct OSFS from url."""
    if url.scheme != "file":
        raise ValueError("Scheme must be == 'file'")
    return OSFS("/")


def make_fs(url):
    """Create PyFilesystem FS for the given url."""
    factories = {"file": make_osfs}
    if url.scheme not in factories:
        raise ValueError("Invalid scheme '%s'" % url.scheme)
    else:
        return factories[url.scheme](url).opendir(url.path)


def fake_data():
    """Create fake ``MetaData`` to make Dash validation happy."""
    return MetaData(id=FAKE_DATA_ID, path="file:///path/to/fake.bam", sample="fake")


def strip_sample(sample):
    """Postprocess sample name"""
    return re.sub(settings.SAMPLE_STRIP_RE, "", sample)


def load_data(url_bam):
    """Load ``MetaData`` from the given ``url_bam``."""
    if url_bam.scheme != "file":
        raise ExcovisException("Can only load file resources at the moment")
    with pysam.AlignmentFile(url_bam.path, "rb") as samfile:
        read_groups = samfile.header.as_dict().get("RG", [])
        if len(read_groups) != 1:
            raise ExcovisException("Must have one read group per BAM file!")
        sample = read_groups[0].get("SM", fs.path.basename(url_bam.path[: -len(".bam")]))
        hash = hashlib.sha256(url_bam.path.encode("utf-8")).hexdigest()
        return MetaData(id=hash, path=url_bam.path, sample=strip_sample(sample))
