# Performance pass

This document describes a set of low-risk performance optimizations added to the
`cellxm` pipeline. All three changes are **transparent by default**: with the
shipped `config.toml`, the output (`cells.gpkg`, `zones.gpkg`, `stations.gpkg`)
is identical to before. The speed levers are opt-in.

## Summary of changes

| # | Change | File(s) | Default impact | Speed lever |
|---|--------|---------|----------------|-------------|
| 1 | Cache the Geofabrik index on disk | `cellxm/io.py` | Output identical | Automatic |
| 2 | Expose a DEM downsample factor | `cellxm/io.py`, `config.toml`, CLI | Output identical (factor=1.0) | `dem_resample_factor` |
| 3 | Persistent multiprocessing pool | `cellxm/main.py`, `config.toml` | Output identical | `maxtasksperchild` |

---

## 1. Geofabrik index caching

`download_osm_elems` previously fetched and parsed the multi-megabyte
`https://download.geofabrik.de/index-v1.json` on **every run**. In a full
`run.sh` over 150+ countries this is the same file downloaded 150+ times.

The index is now fetched through `get_geofabrik_index(...)`, which caches the raw
JSON under `./.cache/osm/geofabrik-index-v1.json` and only hits the network on a
cache miss or once the cache is older than `ttl_hours` (default 720h / 30 days).

Robustness details:
- The cache is written atomically (`*.json.tmp` then `replace`), so an
  interrupted run cannot leave a half-written cache.
- If the network is unavailable but a cache file exists (even a stale one), the
  cached copy is used instead of failing.
- A corrupt/unreadable cache transparently triggers a re-download.

Force a refresh at any time by deleting the cache file, or call
`get_geofabrik_index(..., update=True)`.

## 2. DEM downsample factor (`dem_resample_factor`)

The Copernicus GLO-90 DEM is the input to the topographic analysis, which is the
compute-heavy part of the pipeline. The code already contained a (disabled)
resampling path hardcoded to `1/1`. It is now a real config / CLI option.

```toml
# config.toml
dem_resample_factor=1.0  # 1.0 = no downsampling (default, unchanged behaviour)
```

```bash
# downsample the DEM to ~1/4 the points (0.5 per axis) for a large topo speedup
python cellxm/main.py locate GR --dem-resample-factor 0.5 --out-folder ./data/outputs/GR
```

- `1.0` (or any value `>= 1`) is a no-op — identical to the original pipeline.
- `0.5` halves each axis, i.e. ~4x fewer elevation points. This is the single
  biggest CPU reduction for the topographic step, at a modest accuracy cost in
  fine terrain detail. Recommended starting point when you want speed.
- The factor never upscales and always keeps at least one pixel per axis.

The geometry of the downsample is computed by the pure helper
`io.resample_out_shape(count, height, width, factor)` (unit-tested).

> Note: because this changes how finely topography is sampled, it **can** change
> station counts in mountainous cells. Keep it at `1.0` for the canonical,
> reward-relevant output; use `0.5` for fast exploratory or regional runs and
> diff the results before adopting.

## 3. Persistent multiprocessing pool

Previously a `spawn` pool was created and destroyed **inside** the per-batch
loop, repaying worker-startup cost (which re-imports all modules in each worker)
on every batch. The pool is now created **once** and reused across all batches.

To bound memory (the original per-batch pattern was a workaround for occasional
freezing / leaks), workers are recycled with `maxtasksperchild`:

```toml
# config.toml
maxtasksperchild=200  # recycle each worker after this many tiles
```

Lower this value if you see memory growth on very large countries; raise it (or
set it higher) to maximise reuse on smaller runs.

---

## A note on "S3"

The pipeline downloads elevation tiles from OpenTopography's **public** S3
bucket (`opentopography.s3.sdsc.edu`) anonymously (`no_sign_request=True`). This
is just the data CDN for the Copernicus GLO-90 / SRTM / NASADEM rasters — it is
**not** WeatherXM infrastructure and needs no AWS credentials. Downloaded tiles
are cached locally under `./.cache`, so subsequent runs reuse them. The whole
pipeline runs locally; "needs S3" simply means "needs to download the global DEM".

The dominant cost of a full **world** run is data volume (downloading global DEM
tiles + emitting tens of millions of features), not CPU. For per-country or
per-region runs the CPU optimizations above (especially `dem_resample_factor`)
make a large difference; for a full planet run, parallelise across regions and
pre-warm the `.cache`.

---

## Testing

Unit tests for the pure logic of this pass live in
`cellxm/tests/test_perf.py`. They stub the heavy GDAL/geopandas/rasterio stack so
they run without a full GIS install:

```bash
# as a standalone script (no pytest needed)
python cellxm/tests/test_perf.py

# or under pytest
pytest cellxm/tests/test_perf.py
```

Coverage:
- `resample_out_shape`: no-op cases, half/quarter sizing, min-1-pixel flooring.
- `get_geofabrik_index`: downloads once then reuses, `update=True` refetch,
  TTL expiry, network-outage fallback to cache, corrupt-cache recovery.
- `overwrite_config`: the `dem_resample_factor` plumbing.
