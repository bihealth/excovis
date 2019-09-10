"""Registeration of Dash app callbacks.

Each ``register_*`` function sets up callbacks for one aspect of the app.  This module does not
any component building itself.  Instead, this is done in the module ``.ui``.
"""

import base64
from io import BytesIO

# from plotly.tools import mpl_to_plotly
import dash
import dash_html_components as html
from logzero import logger

from . import plot


def fig_to_uri(in_fig, **save_args):
    out_img = BytesIO()
    in_fig.savefig(out_img, format="png", **save_args)
    out_img.seek(0)  # rewind file
    encoded = base64.b64encode(out_img.read()).decode("ascii").replace("\n", "")
    return "data:image/png;base64,{}".format(encoded)


def register_plot(app):
    """Register the display of the coverage plot."""

    @app.callback(
        dash.dependencies.Output("page-plot", "children"),
        [dash.dependencies.Input("input_%s" % s, "value") for s in ("padding", "gene", "samples")],
    )
    def render_plot(padding, gene, samples):
        logger.info("padding=%s, gene=%s, samples=%s", padding, gene, samples)
        if not gene or not samples:
            return html.Div(
                "After selecting gene and sample(s), the coverage plot will appear here.",
                className="text-center",
            )
        else:
            mpl_fig = plot.render_plot(padding, gene, samples)
            return html.Img(id="cov-plot", src=fig_to_uri(mpl_fig))
