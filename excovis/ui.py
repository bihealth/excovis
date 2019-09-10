"""Definitions of Dash layout."""

import os.path

import dash_bootstrap_components as dbc
import dash_core_components as dcc
import dash_html_components as html

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
        {"label": symbol, "value": symbol} for symbol in sorted(genes.load_transcripts().keys())
    ]
    samples_options = [
        {"label": record.sample, "value": record.id} for record in store.load_all_data()
    ]
    return html.Div(
        children=[
            dbc.Label("Select Gene", html_for="id_input_gene"),
            dcc.Dropdown(id="id_input_gene", options=genes_options),
            dbc.Label("Select Sample(s)", html_for="id_input_samples", className="pt-3"),
            dcc.Dropdown(id="id_input_sample", options=samples_options, multi=True),
        ],
        id="menu",
    )


def render_main_content():
    """Render page main content"""
    return html.Div(
        children=[
            dbc.Row(
                [
                    dbc.Col(
                        # content will be rendered in this element
                        children=render_form(),
                        className="col-2",
                    ),
                    dbc.Col(
                        # content will be rendered in this element
                        children=[html.Div(id="page-content")],
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


def render_dataset_page(data):
    """Render the page main content for dataset visualization."""
    return dbc.Tabs(
        children=[
            # dbc.Tab(
            #     html.Div(
            #         dbc.Row(dbc.Col(dcc.Markdown(data.metadata.readme))), className="mx-2 mt-2"
            #     ),
            #     label="About",
            #     tab_id="tab-about",
            # ),
            # dbc.Tab(
            #     html.Div(cells.render(data), className="mx-2 mt-2"),
            #     label="Cell Annotation",
            #     tab_id="tab-cells",
            # ),
            # dbc.Tab(
            #     html.Div(genes.render(data), className="mx-2 mt-2"),
            #     label="Gene Expression",
            #     tab_id="tab-genes",
            # ),
        ],
        active_tab="tab-cells",
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
