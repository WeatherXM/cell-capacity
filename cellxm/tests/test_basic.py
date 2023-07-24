import tempfile
import pathlib

from typer.testing import CliRunner
from cellxm.main import app

runner = CliRunner()


def test_app():
    with tempfile.TemporaryDirectory() as tdir:
        tdir = pathlib.Path(tdir)
        result = runner.invoke(
            app, ["locate", "AD", "--h3-resolution", "5", "--out-folder", f"{tdir}"]
        )
        assert result.exit_code == 0
        assert "2" in result.stdout
        assert "10" in result.stdout
