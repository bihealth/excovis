"""Registeration of Dash app callbacks.

Each ``register_*`` function sets up callbacks for one aspect of the app.  This module does not
any component building itself.  Instead, this is done in the module ``.ui``.
"""

import base64
from io import BytesIO

import dash
import dash_html_components as html
import dash_table
from intervaltree import Interval, IntervalTree
from logzero import logger
import pandas as pd
from natsort import natsorted

from . import plot, genes, store


def fig_to_uri(in_fig, **save_args):
    out_img = BytesIO()
    in_fig.savefig(out_img, format="png", **save_args)
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)


def register_transcript_select(app):
    """Register updating of transcript selection."""

    @app.callback(
        [
            dash.dependencies.Output("input_transcript", "options"),
            dash.dependencies.Output("input_transcript", "value"),
        ],
        [dash.dependencies.Input("input_gene", "value")],
    )
    def transcript_options(gene_symbol):
        if not gene_symbol:
            return [], []
        else:
            options = [
                {"label": tx.tx_accession, "value": tx.tx_accession}
                for tx in natsorted(
                    filter(
                        lambda tx: tx.gene_symbol == gene_symbol, genes.load_transcripts().values()
                    ),
                    key=lambda tx: tx.tx_accession,
                )
            ]
            return options, options[0]["value"]


def register_plot(app):
    """Register the display of the coverage plot."""

    @app.callback(
        dash.dependencies.Output("page-plot", "children"),
        [
            dash.dependencies.Input("input_%s" % s, "value")
            for s in ("padding", "ymax", "transcript", "samples")
        ],
    )
    def render_plot(padding, ymax, gene, samples):
        if not gene or not samples:
            return html.Div(
                "After selecting gene and sample(s), the coverage plot will appear here.",
                className="text-center",
            )
        else:
            mpl_fig = plot.render_plot(padding, ymax, gene, samples)
            return html.Img(id="cov-plot", src=fig_to_uri(mpl_fig))


def register_table(app):
    """Register the display of the coverage table."""

    @app.callback(
        dash.dependencies.Output("page-table", "children"),
        [
            dash.dependencies.Input("input_%s" % s, "value")
            for s in ("transcript", "samples", "aggregation")
        ],
    )
    def render_plot(tx_accession, samples, aggregation):
        if not tx_accession or not samples:
            return []
        else:
            transcript = genes.load_transcripts()[tx_accession]
            tree = IntervalTree([Interval(exon.begin, exon.end) for exon in transcript.exons])
            coverage_df = store.load_coverage_df(0, tx_accession, samples)
            # Filter for on-target coverage only
            coverage_df = coverage_df[
                coverage_df.apply(lambda x: tree.overlaps_point(x["pos"] - 1), axis=1)
            ].drop(columns=["chrom", "pos"])
            # Make data tidy
            coverage_df = pd.melt(
                frame=coverage_df, id_vars=["exon_no"], var_name="sample", value_name="coverage"
            )
            orig_coverage_df = coverage_df.copy(deep=True)
            # Aggregate per sample and exon.
            grouped_df = coverage_df.groupby(["exon_no", "sample"])
            if aggregation == "min":
                coverage_df = grouped_df.min()
            elif aggregation == "max":
                coverage_df = grouped_df.max()
            elif aggregation == "median":
                coverage_df = grouped_df.median()
            else:  # if aggregation == "mean":
                coverage_df = grouped_df.mean()
            # Make each sample its own column again
            coverage_df = coverage_df.reset_index()
            coverage_df["coverage"] = coverage_df["coverage"].apply(lambda x: round(x, 1))
            coverage_df = coverage_df.pivot(
                index="exon_no", columns="sample", values="coverage"
            ).reset_index()
            coverage_df.insert(0, "feature", "exon")
            # Aggregate per sample.
            grouped_df2 = orig_coverage_df.groupby(["sample"])
            if aggregation == "min":
                coverage_df2 = grouped_df2.min()
            elif aggregation == "max":
                coverage_df2 = grouped_df2.max()
            elif aggregation == "median":
                coverage_df2 = grouped_df2.median()
            else:  # if aggregation == "mean":
                coverage_df2 = grouped_df2.mean()
            # Make each sample its own column again
            coverage_df2 = coverage_df2.reset_index()
            coverage_df2["coverage"] = coverage_df2["coverage"].apply(lambda x: round(x, 1))
            coverage_df2 = coverage_df2.pivot(columns="sample", values="coverage")
            coverage_df2.insert(0, "exon_no", None)
            coverage_df2.insert(0, "feature", "transcript")
            return [
                html.H3("Coverage Table (%s)" % aggregation),
                dash_table.DataTable(
                    columns=[{"name": i, "id": i} for i in coverage_df.columns],
                    data=pd.concat([coverage_df2, coverage_df]).to_dict("records"),
                    sort_action="native",
                ),
            ]
