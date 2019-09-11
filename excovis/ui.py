"""Definitions of Dash layout."""

import os.path

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html
from natsort import natsorted

from . import settings, store, genes
from .__init__ import __version__


def render_dataset(identifier):
    """Display the dataset."""
    data = store.load_data(identifier)
    if data:
        return render_dataset_page(data)
    else:
        return ""


def render_navbar():
    """Render the site navbar"""
    return dbc.Navbar(
        dbc.Container(
            children=[
                # Use row and col to control vertical alignment of logo / brand
                dbc.NavbarBrand(
                    [html.I(className="fas fa-chart-area mr-1"), settings.HOME_BRAND],
                    id="page-brand",
                )
            ]
        ),
        # className="mb-5",
        color="primary",
        dark=True,
        id="page-navbar",
    )


def render_form():
    """Render form for selecting genes and samples."""
    genes_options = [
        {"label": gene_symbol, "value": gene_symbol}
        for gene_symbol in natsorted(
            set((tx.gene_symbol for tx in genes.load_transcripts().values()))
        )
    ]
    samples_options = [
        {"label": record.sample, "value": record.id} for record in store.load_all_data()
    ]
    return html.Div(
        children=[
            dbc.Label("Exon Padding", html_for="input_padding"),
            dcc.Slider(
                id="input_padding",
                value=settings.DEFAULT_EXON_PADDING,
                min=0,
                max=settings.MAX_EXON_PADDING,
                marks={i: "%dbp" % i for i in range(0, settings.MAX_EXON_PADDING + 1, 25)},
            ),
            dbc.Label("Max coverage", html_for="input_ymax", className="pt-3 mt-3"),
            dcc.Slider(
                id="input_ymax",
                value=settings.DEFAULT_MAX_COVERAGE,
                min=0,
                max=settings.MAX_MAX_COVERAGE,
                marks={i: "%dx" % i for i in range(0, settings.MAX_MAX_COVERAGE + 1, 50)},
                className="pb-3",
            ),
            html.Hr(),
            dbc.Label("Select Gene", html_for="input_gene"),
            dcc.Dropdown(id="input_gene", options=genes_options),
            dbc.Label("Select Transcript", html_for="input_transcript", className="pt-3"),
            dcc.Loading(children=[dcc.Dropdown(id="input_transcript")]),
            html.Hr(),
            dbc.Label("Select Sample(s)", html_for="input_samples"),
            dcc.Dropdown(id="input_samples", options=samples_options, multi=True),
            html.Hr(),
            dbc.Label("Coverage Aggregation", html_for="input_aggregation"),
            dcc.Dropdown(
                id="input_aggregation",
                options=[{"label": agg, "value": agg} for agg in ("min", "max", "median", "mean")],
                value="mean",
            ),
        ],
        id="menu",
    )


def render_main_content():
    """Render page main content"""
    return html.Div(
        children=[
            dbc.Row(
                [
                    dbc.Col(children=render_form(), className="col-2"),
                    dbc.Col(
                        # content will be rendered in this element
                        children=[
                            dcc.Loading(children=[html.Div(id="page-plot")]),
                            dcc.Loading(children=[html.Div(id="page-table")]),
                        ],
                        className="col-10",
                    ),
                ]
            )
        ],
        className="container pt-3",
    )


def render_footer():
    """Render page footer"""
    return html.Footer(
        html.Div(
            children=[
                html.Div(
                    children=[
                        html.Div(
                            children=[
                                html.Span(
                                    "ExCoVis v%s by BIH CUBI" % __version__, className="text-muted"
                                )
                            ],
                            className="col-6",
                        ),
                        html.Div(
                            children=[
                                html.A(
                                    children=[
                                        html.I(className="fas fa-globe-europe mr-1"),
                                        "CUBI Homepage",
                                    ],
                                    href="https://www.cubi.bihealth.org",
                                    className="text-muted mr-3",
                                ),
                                html.A(
                                    children=[
                                        html.I(className="fab fa-github mr-1"),
                                        "GitHub Project",
                                    ],
                                    href="https://github.com/bihealth/excovis",
                                    className="text-muted",
                                ),
                            ],
                            className="col-6 text-right",
                        ),
                    ],
                    className="row",
                )
            ],
            className="container",
        ),
        className="footer",
    )


def build_layout():
    """Build the overall Dash app layout"""
    return html.Div(
        children=[
            # Represents the URL bar, doesn't render anything.
            dcc.Location(id="url", refresh=False),
            # Navbar, content, footer.
            render_navbar(),
            render_main_content(),
            render_footer(),
        ],
        id="_dash-app-content",  # make pytest-dash happy
    )
