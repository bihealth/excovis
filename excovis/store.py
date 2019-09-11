"""Access to all data outside of the session.

This module uses the functions from ``data`` to load the data from the appropriate location
and memoizes the loaded data in the Flask cache.  This module is aware of the Flask cache and
thus also of the Dash app.

Note well that the behaviour of iRODS is special because of how ticket access is implemented.
We have to get direct access to the collection for which we have a ticket.  Then we need to
perform a linear search for the collection's data objects and only THEN can we open them.
"""

from itertools import chain

from logzero import logger
import pandas as pd
from intervaltree import Interval, IntervalTree
import pysam

from . import data, genes, settings
from .exceptions import ExcovisException
from .cache import cache


@cache.memoize()
def load_all_data():
    """Load all meta data information from ``settings.DATA_SOURCES``.

    A data source can either be a URL to a file ending on ``.bam`` or a directory that contains ``.bam`` files.
    """
    result = []

    if settings.FAKE_DATA:
        result.append(data.fake_data())

    for url in settings.DATA_SOURCES:
        if url.scheme in data.PYFS_SCHEMES:
            if url.path.endswith(".bam"):  # one file
                result.append(data.load_data(url))
            else:
                curr_fs = data.make_fs(url)
                for match in curr_fs.glob("**/*.bam"):
                    x = url._replace(path=url.path + match.path)
                    result.append(data.load_data(x))
    return result


@cache.memoize()
def load_data(id):
    for data in load_all_data():
        if data.id == id:
            return data
    raise ExcovisException("Unknown dataset %d" % id)


def _load_fake_coverage(sample, chrom, tree):
    def padded_range(a, b, padding):
        return range(a - padding, b + padding)

    def fn(lst):
        return list(sorted(set(chain(*lst))))

    positions = fn([padded_range(itv.begin, itv.end, settings.MAX_EXON_PADDING) for itv in tree])
    n = len(positions)
    return pd.DataFrame(
        data=[
            {"chrom": chrom, "pos": pos, sample: int(50.0 * i / n)}
            for i, pos in enumerate(positions)
        ],
        columns=["chrom", "pos", sample],
    )


@cache.memoize()
def load_coverage(sample_id, chrom, tree, transcript):
    """Load coverage for all positions in ``tree`` from ``chrom``."""
    if sample_id == data.FAKE_DATA_ID:  # short-circuit for fake data
        return _load_fake_coverage(sample_id, chrom, tree)

    datasets = load_all_data()
    for dataset in datasets:
        if dataset.id == sample_id:
            break
    else:
        logger.info("Could not locate sample %s in %s", sample_id, [ds.id for ds in datasets])
        raise ExcovisException("Unknown sample %s" % sample_id)

    logger.info("dataset = %s", dataset)

    pad = settings.MAX_EXON_PADDING
    rows = []

    with pysam.AlignmentFile(dataset.path, "rb") as samfile:
        for i, itv in enumerate(sorted(tree, key=lambda exon: exon.begin)):
            if transcript.strand == "+":
                exon_no = i + 1
            else:
                exon_no = len(transcript.exons) - i
            seen = set()
            for align_col in samfile.pileup(chrom, itv.begin - pad, itv.end + pad):
                pos = align_col.reference_pos
                if pos not in seen and itv.begin - pad <= pos < itv.end + pad:
                    seen.add(pos)
                    rows.append(
                        {
                            "chrom": chrom,
                            "pos": pos + 1,
                            "exon_no": exon_no,
                            dataset.sample: align_col.get_num_aligned(),
                        }
                    )
            for pos in range(itv.begin - pad, itv.end + pad):
                if pos not in seen:
                    rows.append(
                        {"chrom": chrom, "pos": pos + 1, "exon_no": exon_no, dataset.sample: 0}
                    )
    result = pd.DataFrame(data=rows, columns=["chrom", "pos", "exon_no", dataset.sample])
    result.sort_values("pos", inplace=True)
    return result


@cache.memoize()
def load_coverage_df(exon_padding, tx_accession, samples):
    transcript = genes.load_transcripts()[tx_accession]
    tree = IntervalTree([Interval(exon.begin, exon.end) for exon in transcript.exons])
    ds = [load_coverage(sample, transcript.chrom, tree, transcript) for sample in samples]
    df_coverage = pd.concat(
        [ds[0]["chrom"], ds[0]["pos"], ds[0]["exon_no"]] + [d.iloc[:, 3] for d in ds],
        axis="columns",
    )
    df_coverage.sort_values("pos", inplace=True)
    return df_coverage
