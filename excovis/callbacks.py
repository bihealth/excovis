"""Registeration of Dash app callbacks.

Each ``register_*`` function sets up callbacks for one aspect of the app.  This module does not
any component building itself.  Instead, this is done in the module ``.ui``.
"""

import base64
import os.path
import uuid

import dash
import dash_html_components as html
from logzero import logger
from werkzeug.utils import secure_filename

from . import ui, settings, store


# def register_page_content(app):
#     """Register the display of the page content with the app."""
#
#     @app.callback(
#         dash.dependencies.Output("page-content", "children"),
#         [dash.dependencies.Input("url", "pathname")],
#     )
#     def render_page_content(pathname):
#         view, kwargs = get_route(pathname)
#         if view == "home":
#             return ui.render_home()
#         elif view == "viz":
#             return ui.render_dataset(kwargs.get("dataset"))
#         else:
#             return ui.render_not_found()
#
#
# def register_page_brand(app):
#     """Register the display of the page brand with the app."""
#
#     @app.callback(
#         dash.dependencies.Output("page-brand", "children"),
#         [dash.dependencies.Input("url", "pathname")],
#     )
#     def render_page_brand(pathname):
#         view, kwargs = get_route(pathname)
#         if view == "home":
#             return [html.I(className="fas fa-home mr-1"), settings.HOME_BRAND]
#         elif view == "viz":
#             metadata = store.load_data(kwargs.get("dataset"))
#             if metadata:
#                 return [html.I(className="fas fa-file-alt mr-1"), metadata.sample]
#             else:
#                 return [html.I(className="fas fa-file-alt mr-1"), "Not Found"]
#         else:
#             return [html.I(className="fas fa-file-alt mr-1"), "Not Found"]
