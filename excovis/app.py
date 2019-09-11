"""Setup the ExCoVis Dash application.

When importing this module, the app is built and configured.  Thus, it is important that before
this module is imported, the values in ``.settings`` must already have been setup.
"""

import os

import dash
import flask


from . import cache, callbacks, settings
from .__init__ import __version__
from .ui import build_layout

#: Path to assets.
ASSETS_FOLDER = os.path.join(os.path.dirname(__file__), "assets")

#: The Flask application to use.
app_flask = flask.Flask(__name__)

# Setup temporary upload folder
app_flask.config["UPLOAD_FOLDER"] = settings.TEMP_DIR
# Setup maximal file upload size
app_flask.config["MAX_CONTENT_LENGTH"] = settings.MAX_UPLOAD_SIZE
# Setup URL prefix for Flask.
app_flask.config["APPLICATION_ROOT"] = "%s/" % settings.PUBLIC_URL_PREFIX

#: The Dash application to run.
app = dash.Dash(
    __name__,
    # Use our specific Flask app
    server=app_flask,
    # All files from "assets/" will be served as "/assets/*"
    assets_folder=ASSETS_FOLDER,
    # The visualization will be served below "/dash"
    routes_pathname_prefix="/dash/",
    requests_pathname_prefix="%s/dash/" % settings.PUBLIC_URL_PREFIX,
)

# Setup the cache.
cache.setup_cache(app)

# Set app title
app.title = "ExCoVis v%s" % __version__

# Serve assets locally
app.css.config.serve_locally = True
app.scripts.config.serve_locally = True

# Setup the application's main layout.
app.layout = build_layout()

# Register the callbacks with the app.
#
# Callbacks for the coverage plot.
callbacks.register_plot(app)
callbacks.register_table(app)
callbacks.register_transcript_select(app)

# Add redirection for root.
@app_flask.route("/")
def redirect_root():
    return flask.redirect("%s/dash/" % settings.PUBLIC_URL_PREFIX)
