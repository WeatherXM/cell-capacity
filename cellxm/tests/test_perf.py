"""Unit tests for the performance-pass changes.

These tests cover the *pure*, deterministic logic added in the performance pass:
  * io.resample_out_shape   - DEM downsample geometry
  * io.get_geofabrik_index  - on-disk caching of the Geofabrik index
  * utils.overwrite_config  - config plumbing of the new lever

The heavy geospatial C-extension dependencies (GDAL/rasterio/geopandas/...) are
stubbed out so the logic can be exercised without a full GIS install. This keeps
the tests hermetic and fast; the real pipeline behaviour is unchanged because
these functions are imported from the actual source modules.
"""

import sys
import types
import json
import time
import tempfile
import pathlib


# --------------------------------------------------------------------------- #
# Stub the heavy / unavailable third-party modules before importing cellxm.*   #
# --------------------------------------------------------------------------- #
def _stub(name, **attrs):
    if name in sys.modules:
        mod = sys.modules[name]
    else:
        mod = types.ModuleType(name)
        sys.modules[name] = mod
    for k, v in attrs.items():
        setattr(mod, k, v)
    return mod


class _FakeResponse:
    def __init__(self, status_code=200, payload=None):
        self.status_code = status_code
        self._payload = payload if payload is not None else {"type": "FeatureCollection"}

    def json(self):
        return self._payload


class _FakeRequests:
    """Minimal stand-in for the `requests` module with a call counter."""

    def __init__(self):
        self.calls = 0
        self.next_status = 200
        self.payload = {"type": "FeatureCollection", "tag": "from-network"}

    def get(self, url, *a, **kw):
        self.calls += 1
        return _FakeResponse(self.next_status, self.payload)


_fake_requests = _FakeRequests()
_stub("requests").get = _fake_requests.get
_stub("requests.adapters", HTTPAdapter=object, Retry=object)

# geo / raster stack -> empty stubs (only needed so the import line succeeds)
for _m in [
    "rasterio", "rasterio.features", "rasterio.mask", "rioxarray", "shapely",
]:
    _stub(_m)
_stub("geopandas", GeoDataFrame=object, GeoSeries=object)
_stub("shapely.geometry", Polygon=object, MultiPolygon=object)
_stub("shapely.validation", make_valid=lambda x: x)
_stub("rasterio.enums", Resampling=types.SimpleNamespace(bilinear=1))
_stub("rasterio.merge", merge=lambda *a, **k: (None, None))
_stub("osgeo", gdal=types.SimpleNamespace())
_stub("cloudpathlib", S3Client=object)
# cellxm.ogr2ogr pulls in osgeo at import time -> stub the whole submodule
_stub("cellxm.ogr2ogr", main=lambda *a, **k: 0)

# Make `cellxm` importable as a namespace pointing at the real package dir
_pkg_dir = pathlib.Path(__file__).resolve().parents[1]   # .../cellxm
_repo_root = _pkg_dir.parent
if str(_repo_root) not in sys.path:
    sys.path.insert(0, str(_repo_root))

from cellxm.io import resample_out_shape, get_geofabrik_index, GEOFABRIK_INDEX_URL  # noqa: E402
from cellxm.utils import overwrite_config  # noqa: E402


# --------------------------------------------------------------------------- #
# resample_out_shape                                                           #
# --------------------------------------------------------------------------- #
def test_resample_noop_none():
    assert resample_out_shape(1, 100, 80, None) == (1, 100, 80)

def test_resample_noop_factor_one():
    assert resample_out_shape(1, 100, 80, 1.0) == (1, 100, 80)

def test_resample_noop_factor_gt_one():
    # values >= 1 must never upscale
    assert resample_out_shape(2, 100, 80, 1.5) == (2, 100, 80)

def test_resample_half():
    assert resample_out_shape(1, 100, 80, 0.5) == (1, 50, 40)

def test_resample_quarter_area_is_half_each_axis():
    c, h, w = resample_out_shape(1, 360, 360, 0.5)
    # 0.5 per-axis  ->  ~1/4 the points  (the headline speed win)
    assert h * w == (360 * 360) // 4

def test_resample_floor_to_at_least_one_pixel():
    # tiny rasters / aggressive factors must not collapse to 0 px
    assert resample_out_shape(1, 1, 1, 0.5) == (1, 1, 1)
    assert resample_out_shape(1, 10, 10, 0.01) == (1, 1, 1)


# --------------------------------------------------------------------------- #
# get_geofabrik_index caching                                                 #
# --------------------------------------------------------------------------- #
def _reset_net():
    _fake_requests.calls = 0
    _fake_requests.next_status = 200

def test_index_cache_downloads_once_then_reuses():
    _reset_net()
    with tempfile.TemporaryDirectory() as d:
        a = get_geofabrik_index(d)
        b = get_geofabrik_index(d)          # should hit the on-disk cache
        assert _fake_requests.calls == 1, "second call must not hit the network"
        assert a == b
        assert (pathlib.Path(d) / "geofabrik-index-v1.json").exists()

def test_index_cache_update_forces_refetch():
    _reset_net()
    with tempfile.TemporaryDirectory() as d:
        get_geofabrik_index(d)
        get_geofabrik_index(d, update=True)
        assert _fake_requests.calls == 2

def test_index_cache_ttl_expiry_refetches():
    _reset_net()
    with tempfile.TemporaryDirectory() as d:
        get_geofabrik_index(d)
        # ttl 0 => existing cache is always considered stale
        get_geofabrik_index(d, ttl_hours=0)
        assert _fake_requests.calls == 2

def test_index_cache_used_when_network_down():
    _reset_net()
    with tempfile.TemporaryDirectory() as d:
        get_geofabrik_index(d)              # populate cache (call 1)
        _fake_requests.next_status = 500    # simulate outage
        out = get_geofabrik_index(d, ttl_hours=0)  # forced refetch fails -> fallback
        assert out is not None
        assert _fake_requests.calls == 2    # it tried, failed, used cache

def test_index_cache_recovers_from_corrupt_file():
    _reset_net()
    with tempfile.TemporaryDirectory() as d:
        get_geofabrik_index(d)              # call 1
        (pathlib.Path(d) / "geofabrik-index-v1.json").write_text("{ not json")
        get_geofabrik_index(d)              # corrupt -> redownload (call 2)
        assert _fake_requests.calls == 2


# --------------------------------------------------------------------------- #
# overwrite_config plumbing                                                    #
# --------------------------------------------------------------------------- #
def _base_args(cfg, **over):
    # positional signature of overwrite_config minus the leading config dict
    defaults = dict(
        h3_resolution=None, elev_diff_thresh=None, cell_buffer=None,
        min_area_urban=None, min_area_green=None, largest_urban=None,
        largest_green=None, min_aspect_zone_area_cell_perc=None,
        include_buildings=None, outfolder=None,
    )
    defaults.update(over)
    return overwrite_config(
        cfg,
        defaults["h3_resolution"], defaults["elev_diff_thresh"],
        defaults["cell_buffer"], defaults["min_area_urban"],
        defaults["min_area_green"], defaults["largest_urban"],
        defaults["largest_green"], defaults["min_aspect_zone_area_cell_perc"],
        defaults["include_buildings"], defaults["outfolder"],
        dem_resample_factor=defaults.get("dem_resample_factor"),
    )

def test_config_resample_override_applied():
    cfg = _base_args({}, dem_resample_factor=0.5)
    assert cfg["dem_resample_factor"] == 0.5

def test_config_resample_absent_when_not_set():
    cfg = _base_args({"dem_resample_factor": 1.0})
    # passing None must not clobber the existing config value
    assert cfg["dem_resample_factor"] == 1.0


# --------------------------------------------------------------------------- #
# Standalone runner (no pytest dependency required)                           #
# --------------------------------------------------------------------------- #
if __name__ == "__main__":
    fns = [v for k, v in sorted(globals().items())
           if k.startswith("test_") and callable(v)]
    passed = 0
    for fn in fns:
        try:
            fn()
            print(f"PASS  {fn.__name__}")
            passed += 1
        except Exception as e:
            print(f"FAIL  {fn.__name__}: {e!r}")
            raise
    print(f"\n{passed}/{len(fns)} tests passed")
