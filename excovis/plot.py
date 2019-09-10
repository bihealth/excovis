"""Code for plotting."""

from itertools import chain
import os
import sys

import pandas as pd
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.path as pltpath
import matplotlib.patches as patches

from excovis.exceptions import ExcovisException
from . import data, genes, settings

from logzero import logger
from intervaltree import Interval, IntervalTree
import pysam
import numpy as np


MARGIN_X = 100
PADDING = 100  # bp
MARGIN_Y_TX = 2
EXON_UTR_HEIGHT = 5
EXON_CDS_HEIGHT = 10
EXON_ROW_HEIGHT = 20
TEXT_MARGIN = 1

#: Maximal coverage to show.
MAX_COV = 50
#: Minimal coverage for level "warning" (smaller triggers error).
MIN_WARN = 10
#: Minimal coverage for level "OK" (smaller triggers warning).
MIN_OK = 20

#: Figure height
FIGSIZE_H = 12
#: Figure widget
FIGSIZE_V = 2.5


def pos_filtered(df, chrom, begin, end):
    return df[(df.chrom == chrom) & (df.pos >= begin - 1) & (df.pos < end)]


def val_to_qual(x):
    """Convert the given quality value to a quality level."""
    if x < MIN_WARN:
        return 0
    elif x < MIN_OK:
        return 1
    else:
        return 2


def build_qualmap(quals):
    """Build a matplotlib color map for the given qualities."""
    lst = []
    if 0 in quals:  # FAIL
        lst.append((1.0, 0.39, 0.34, 1.0))
    if 1 in quals:  # warn
        lst.append((0.98, 0.81, 0.33, 1.0))
    if 2 in quals:  # OK
        lst.append((0.38, 1.0, 0.46, 1.0))
    return mpl.colors.ListedColormap(lst)


def _plot_for_gene(
    transcripts, gene_symbol, df_covs, sep_vlines=[], projected=False, suptitle=None
):
    if transcripts.empty:
        tx_chrom = "1"
        tx_strand = "+"
    else:
        tx_chrom = transcripts.chrom.iloc[0]
        tx_strand = transcripts.strand.iloc[0]

    # Compute start and end genome position.
    pos_begin = transcripts.begin.min() - PADDING - MARGIN_X
    pos_end = transcripts.end.max() + PADDING + MARGIN_X
    # Extract coverage information from coverage data frame.
    tx_covs = pos_filtered(df_covs, tx_chrom, pos_begin, pos_end)
    # Get sample names.
    samples = list(tx_covs.columns[2:])

    # Initialize figure.
    fig = plt.figure(figsize=(FIGSIZE_H, (len(samples) + 1) * FIGSIZE_V), dpi=75)
    if suptitle:
        fig.suptitle(suptitle)

    # Collect CDS begin/end positions
    num_tx = transcripts.shape[0]
    cds_vlines = set()
    for idx in range(num_tx):
        transcript = transcripts.iloc[idx, :]
        if transcript.cds_begin != transcript.cds_end:
            cds_vlines |= {transcript.cds_begin, transcript.cds_end}

    # Create one coverage sub figure for each sample.
    for sample_idx, sample in enumerate(samples):
        ax = fig.add_subplot(len(samples) + 1, 1, sample_idx + 1)
        logger.info("pos_begin=%d pos_end=%d", pos_begin, pos_end)
        ax.set_xlim(pos_begin, pos_end)
        x = tx_covs.loc[:, sample].apply(val_to_qual).values
        # background image
        im = ax.imshow(
            x.reshape(1, -1),
            cmap=build_qualmap(x),
            aspect="auto",
            origin="lower",
            extent=[tx_covs.pos.min(), tx_covs.pos.max(), 0, MAX_COV],
        )
        # path for masking
        paths = ax.fill_between(
            x=list(tx_covs.pos.values),
            y1=list(tx_covs.loc[:, sample].values),
            facecolor="none",
            edgecolor="none",
        )
        logger.info("%s %s", sample, tx_covs.columns)
        # Make the 'fill' mask and clip the background image with it.
        patch = patches.PathPatch(paths._paths[0], visible=False)
        ax.add_artist(patch)
        im.set_clip_path(patch)
        ax.set(ylabel="depth of coverage", ylim=(0, MAX_COV))
        if projected:
            ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)
        # Show lines with warning and error threshold.
        ax.axhline(y=MIN_OK, lw=1, ls=":", color="black")
        ax.axhline(y=MIN_WARN, lw=1, ls=":", color="black")
        # Display vertical lines indicating CDS start and end
        for x in cds_vlines:
            ax.axvline(x=x, lw=1, ls=":", color="dimgray")
        if sep_vlines:
            # Display vertical lines indicating exon jumps.
            for x in sep_vlines:
                ax.axvline(x=x, lw=1, ls="-", color="dimgray")
        # Show sample name in top right
        ax.text(
            0.99,
            1 - 0.01 * FIGSIZE_H / FIGSIZE_V,
            "sample: %s" % sample,
            horizontalalignment="right",
            verticalalignment="top",
            transform=ax.transAxes,
        )

    # Create bottom plot with the transcripts.
    ax = fig.add_subplot(len(samples) + 1, 1, len(samples) + 1)
    ax.set_yticks([])
    ax.set_xlim(pos_begin, pos_end)
    ax.set_ylim(-num_tx * EXON_ROW_HEIGHT - MARGIN_Y_TX, MARGIN_Y_TX)
    ax.tick_params(axis="y", which="both", bottom=False, top=False, labelbottom=False)
    if sep_vlines:
        # Display vertical lines indicating exon jumps.
        for x in sep_vlines:
            ax.axvline(x=x, lw=1, ls="-", color="dimgray")
    if projected:
        ax.tick_params(axis="x", which="both", bottom=False, top=False, labelbottom=False)

    for idx in range(num_tx):
        transcript = transcripts.iloc[idx, :]
        ax.text(  # transcript name
            pos_begin + (pos_end - pos_begin) / 2,
            -idx * EXON_ROW_HEIGHT - EXON_CDS_HEIGHT - TEXT_MARGIN,
            transcript.tx_accession,
            horizontalalignment="center",
            verticalalignment="top",
            size="larger",
        )
        if not projected:
            ax.add_patch(  # dashed line connecting exons
                patches.PathPatch(
                    pltpath.Path(
                        (
                            (transcript.begin, -idx * EXON_ROW_HEIGHT - 0.5 * EXON_CDS_HEIGHT),
                            (transcript.end, -idx * EXON_ROW_HEIGHT - 0.5 * EXON_CDS_HEIGHT),
                        ),
                        (pltpath.Path.MOVETO, pltpath.Path.LINETO),
                    ),
                    ls="--",
                )
            )
        for begin, length in zip(transcript.exon_begins, transcript.exon_lengths):
            ax.add_patch(  # box with UTR height
                patches.Rectangle(
                    (begin, -idx * EXON_ROW_HEIGHT - EXON_CDS_HEIGHT + EXON_UTR_HEIGHT / 2),
                    length,
                    EXON_UTR_HEIGHT,
                    linewidth=1,
                    color="#3a3a3a",
                )
            )
            cds_begin = transcript.cds_begin
            cds_end = transcript.cds_end
            if cds_begin <= (begin + length) and begin <= cds_end:
                off_left = 0 if begin >= cds_begin else cds_begin - begin
                off_right = 0 if cds_end > (begin + length) else (begin + length) - cds_end
                ax.add_patch(  # box with CDS height
                    patches.Rectangle(
                        (begin + off_left, -idx * EXON_ROW_HEIGHT - EXON_CDS_HEIGHT),
                        length - off_left - off_right,
                        EXON_CDS_HEIGHT,
                        linewidth=1,
                        color="#3a3a3a",
                    )
                )
    fig.subplots_adjust(top=0.95)
    return fig


def _plot_for_gene_projected(transcripts, gene_symbol, df_covs, exon_padding):
    # Prepare the projection from chromosome to plotted space (only consider +/- exon_padding bases around the exon).
    positions = set()
    for transcript in transcripts:
        for exon in transcript.exons:
            positions |= set(range(exon.begin - exon_padding, exon.end + exon_padding + 1))
    projection = {g: p for p, g in enumerate(sorted(positions))}
    # Project transcripts.
    proj_transcripts = pd.DataFrame(
        data=[
            {
                "chrom": transcripts[0].chrom,
                "strand": transcripts[0].strand,
                "tx_accession": transcripts[0].tx_accession,
                "gene_symbol": transcripts[0].gene_symbol,
                "begin": tx.tx_begin,
                "end": tx.tx_end,
                "cds_begin": tx.cds_begin,
                "cds_end": tx.cds_end,
                "exon_begins": [exon.begin for exon in tx.exons],
                "exon_lengths": [exon.length() for exon in tx.exons],
            }
            for tx in transcripts
        ],
        columns=[
            "chrom",
            "gene_symbol",
            "tx_accession",
            "strand",
            "begin",
            "end",
            "cds_begin",
            "cds_end",
            "exon_begins",
            "exon_lengths",
        ],
    )

    def fn(x):
        return projection.get(x)

    def fn2(xs):
        return [projection.get(x) for x in xs]

    proj_transcripts.begin = proj_transcripts.loc[:, "begin"].apply(fn)
    proj_transcripts.end = proj_transcripts.loc[:, "end"].apply(fn)
    proj_transcripts.cds_begin = proj_transcripts.loc[:, "cds_begin"].apply(fn)
    proj_transcripts.cds_end = proj_transcripts.loc[:, "cds_end"].apply(fn)
    proj_transcripts.exon_begins = proj_transcripts.loc[:, "exon_begins"].apply(fn2)

    # Project coverage positions (-1 marks null)
    proj_covs = df_covs.copy()
    proj_covs["pos"] = proj_covs.loc[:, "pos"].apply(lambda pos: projection.get(pos, -1))
    proj_covs = proj_covs[proj_covs.pos >= 0]

    # Compute positions of vertical lines indicating a jump in the coordinate system.
    jump_positions = []
    prev = None
    for gpos, ppos in projection.items():
        if prev is not None and gpos != prev + 1:
            jump_positions.append(ppos)
        prev = gpos

    # Plot.
    logger.info("proj_transcripts %s", proj_transcripts.columns)
    logger.info("proj_covs %s", proj_covs.columns)
    return _plot_for_gene(
        proj_transcripts,
        gene_symbol,
        proj_covs,
        sep_vlines=jump_positions,
        projected=True,
        suptitle="Coverage for gene %s with padding %d bp" % (gene_symbol, exon_padding),
    )


def plot_for_gene(transcripts, gene_symbol, df_covs, exon_padding=None):
    if exon_padding is None:
        # Select transcripts that we are interested in.
        return _plot_for_gene(
            transcripts, gene_symbol, df_covs, suptitle="Coverage for gene %s" % gene_symbol
        )
    else:
        return _plot_for_gene_projected(transcripts, gene_symbol, df_covs, exon_padding)


def load_coverage(sample, chrom, tree):
    """Load coverage for all positions in ``tree`` from ``chrom``."""
    if sample == data.FAKE_DATA_ID:

        def padded_range(a, b, padding):
            return range(a - padding, b + padding)

        def fn(lst):
            return list(sorted(set(chain(*lst))))

        positions = fn(
            [padded_range(itv.begin, itv.end, settings.MAX_EXON_PADDING) for itv in tree]
        )
        n = len(positions)
        return pd.DataFrame(
            data=[
                {"chrom": chrom, "pos": pos, sample: int(50.0 * i / n)}
                for i, pos in enumerate(positions)
            ],
            columns=["chrom", "pos", sample],
        )
    else:
        raise ExcovisException("Implement me! %s" % sample)


def render_plot(exon_padding, gene_symbol, samples):
    transcripts = genes.load_transcripts().get(gene_symbol, [])[:1]
    if not transcripts:
        df_coverage = pd.DataFrame(data=[], columns=["chrom", "pos"] + samples)
    else:
        tree = IntervalTree(
            [
                Interval(exon.begin, exon.end)
                for transcript in transcripts
                for exon in transcript.exons
            ]
        )
        chrom = transcripts[0].chrom
        total_length = sum([itv.length() for itv in tree])
        logger.info("chrom = %s, total_length = %d, tree = %s", chrom, total_length, tree)
        ds = [load_coverage(sample, chrom, tree) for sample in samples]
        df_coverage = ds[0]
        for d in ds[1:]:
            df_coverage = df_coverage.join(d, on=["chrom", "pos"])
        print(df_coverage)
    return plot_for_gene(transcripts, gene_symbol, df_coverage, exon_padding)
