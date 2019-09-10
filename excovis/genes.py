"""Helpers for retrieving genes information."""

import gzip
import typing
from urllib.request import urlopen

import attr
from logzero import logger

from .cache import cache

#: URL to ``ncbiRefSeq.txt.gz`` file for GRCh37.
NCBI_REF_SEQ_GRCH37 = "http://hgdownload.cse.ucsc.edu/goldenPath/hg19/database/ncbiRefSeq.txt.gz"


@attr.s(auto_attribs=True)
class Exon:
    """Represent exon in ``Gene``."""

    #: 0-based begin position.
    begin: int
    #: end position
    end: int

    def length(self):
        return self.end - self.begin


@attr.s(auto_attribs=True)
class Transcript:
    """Represent the genome."""

    #: The gene symbol.
    gene_symbol: str
    #: The transcript identifier.
    tx_accession: str
    #: Strand the exon is on.
    strand: str
    #: Chromosome name.
    chrom: str
    #: 0-based begin position of gene.
    tx_begin: int
    #: end position of gene.
    tx_end: int
    #: 0-based begin position of CDS.
    cds_begin: int
    #: end position of CDS
    cds_end: int
    #: Exons
    exons: typing.Tuple[Exon]


@cache.memoize()
def load_transcripts(url=NCBI_REF_SEQ_GRCH37):
    result = {}
    logger.info("Opening URL %s...", url)
    # with urlopen(url) as gzf:
    #     logger.info("Opening .gz file...")
    #     with gzip.GzipFile(fileobj=gzf, mode="r") as inputf:
    with gzip.open("/tmp/ncbiRefSeq.txt.gz") as inputf:
        for line in inputf:
            arr = line.decode("utf-8").strip().split("\t")
            begins = list(map(int, arr[9].split(",")[:-1]))
            ends = list(map(int, arr[10].split(",")[:-1]))
            transcript = Transcript(
                gene_symbol=arr[12],
                tx_accession=arr[1],
                strand=arr[3],
                chrom=arr[2][3:],
                tx_begin=int(arr[4]),
                tx_end=int(arr[5]),
                cds_begin=int(arr[6]),
                cds_end=int(arr[7]),
                exons=tuple(Exon(begin=begin, end=end) for begin, end in zip(begins, ends)),
            )
            result.setdefault(transcript.gene_symbol, []).append(transcript)
    return result
