"""Code for the definition of the data and meta data and loading them from disk.

The types of this modules are used commonly in the excovis app but the actual ``load_`` functions
should only be accessed through the ``.store`` module such that they can be cached by the app.
This module is unaware of the Dash app.
"""

import os
import os.path
import contextlib
from urllib.parse import urlunparse as _urlunparse
import shutil
import ssl
from urllib.parse import parse_qs, urlunparse

import attr
import fs.path
import fs.tools
from fs.osfs import OSFS
from fs.tempfs import TempFS
from logzero import logger
import numpy as np
import pandas as pd
import requests

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
    return MetaData(id="fake.bam", path="/path/to/fake.bam", sample="fake")
