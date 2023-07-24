from typing import Tuple

from pathlib import Path

import tempfile
import requests
from requests.adapters import HTTPAdapter, Retry


import rasterio
import rasterio.features
import rasterio.mask

import numpy as np
import pandas as pd
import geopandas as gpd
import shapely
from shapely.geometry import Polygon, MultiPolygon
from shapely.validation import make_valid

import rioxarray

from datetime import timedelta
from osgeo import gdal
from rasterio.enums import Resampling

from cloudpathlib import S3Client

# from fiona.errors import DriverError
from pathlib import Path

from cellxm import ogr2ogr


def download_osm_elems(
    country_code: str,
    cache_folder="./.cache/osm",
    update: bool = False,
) -> Polygon | MultiPolygon:
    cache_folder = Path(cache_folder)

    cache_folder.mkdir(parents=True, exist_ok=True)

    fppbf = cache_folder / f"{country_code}.osm.pbf"
    fpgpkg = cache_folder / f"{country_code}.gpkg"

    r = requests.get("https://download.geofabrik.de/index-v1.json")
    if r.status_code != 200:
        raise ConnectionError(
            """Failure to retrieve data from geofabrik.de
                Please check if https://download.geofabrik.de/index-v1.json is accessible"""
        )

    df = gpd.GeoDataFrame.from_features(r.json())
    # df = pd.DataFrame([f["properties"] for f in r.json()["features"]])
    try:
        c = "iso3166-1:alpha2"
        df = df[df[c].astype(str).str.contains(country_code)].iloc[0]  # type:ignore
        url = df["urls"]["pbf"]
        boundary = make_valid(df.geometry)
    except IndexError:
        raise ValueError("Please verify the 2-letter country code")

    if (not fpgpkg.exists()) or update:
        s = requests.Session()

        retries = Retry(
            total=5, backoff_factor=10, status_forcelist=[500, 502, 503, 504]
        )

        s.mount("https://", HTTPAdapter(max_retries=retries))
        pbf = s.get(url)

        if pbf.status_code != 200:
            raise ValueError(f"Did not manage to download the pbf for: {country_code}")

        with tempfile.NamedTemporaryFile(suffix=".osm.pbf") as f:
            f.write(pbf.content)

            s = "landuse,place,aeroway,leisure,natural,building"
            q = """ landuse IN ('residential', 'industrial', 'commercial', 'retail') OR
                    place IN ('neighbourhood') OR aeroway IN ('aerodrome') OR
                    leisure IN ('park') OR natural IN ('water')
                    OR building IN ('yes', 'residential', 'apartments', 'house', 'detached', 'semi-detached')  """

            # "-simplify",
            # "0.0001",

            r = ogr2ogr.main(
                [
                    "ogr2ogr",
                    "-select",
                    s,
                    "-where",
                    q,
                    "-f",
                    "GPKG",
                    f"{fpgpkg}",
                    f.name,
                    "multipolygons",
                ]
            )

            if r == 0:  # successful extract from osm
                osm = gpd.read_file(fpgpkg, engine="pyogrio")
                fpgpkg.unlink()
                osm.to_file(  # type:ignore
                    f"{fpgpkg}", driver="GPKG", engine="pyogrio", mode="w"
                )
            else:
                raise ValueError(
                    f"Conversion of .osm.pbf to .gpkg failed. Delete the {fppbf} file and retry"
                )
    return boundary  # type:ignore


def download_dem_data(bbox: Tuple[int, int, int, int], datasource: str = "COP90"):
    client = S3Client(
        no_sign_request=True,
        endpoint_url="https://opentopography.s3.sdsc.edu",
        local_cache_dir=".cache",
    )

    minx, miny, maxx, maxy = bbox

    files = []

    match datasource:
        case "SRTM_GL3":
            base = "s3://raster/SRTM_GL3/SRTM_GL3_srtm"
            cp = client.CloudPath(base)

            for x in list(range(minx, maxx)):
                for y in list(range(miny, maxy)):
                    absx = abs(x)
                    absy = abs(y)

                    if x >= 0:
                        xd = "E"
                    else:
                        xd = "W"

                    if y >= 0:  # North
                        yd = "N"
                        if 0 <= y < 30:
                            fp = (
                                cp
                                / f"North/North_0_29/{yd}{str(absy).zfill(2)}{xd}{str(absx).zfill(3)}.tif"
                            )
                        else:
                            fp = (
                                cp
                                / f"North/North_30_60/{yd}{str(absy).zfill(2)}{xd}{str(absx).zfill(3)}.tif"
                            )

                    else:
                        yd = "S"
                        fp = (
                            cp
                            / f"South/{yd}{str(absy).zfill(2)}{xd}{str(absx).zfill(3)}.tif"
                        )

                    files.append(fp)

        case "NASADEM":
            base = "s3://raster/NASADEM/NASADEM_be"
            cp = client.CloudPath(base)

            for x in list(range(minx, maxx)):
                for y in list(range(miny, maxy)):
                    absx = abs(x)
                    absy = abs(y)
                    if x >= 0:  # East
                        xd = "e"
                    else:
                        xd = "w"

                    if y >= 0:  # North
                        yd = "n"
                    else:
                        yd = "s"

                    fp = (
                        cp
                        / f"NASADEM_HGT_{yd}{str(absy).zfill(2)}{xd}{str(absx).zfill(3)}.tif"
                    )
                    files.append(fp)
        case "COP90":
            base = "s3://raster/COP90/COP90_hh"
            cp = client.CloudPath(base)

            for x in list(range(minx, maxx)):
                for y in list(range(miny, maxy)):
                    absx = abs(x)
                    absy = abs(y)
                    if x >= 0:  # East
                        xd = "E"
                    else:
                        xd = "W"

                    if y >= 0:  # North
                        yd = "N"
                    else:
                        yd = "S"

                    fp = (
                        cp
                        / f"Copernicus_DSM_COG_30_{yd}{str(absy).zfill(2)}_00_{xd}{str(absx).zfill(3)}_00_DEM.tif"
                    )
                    files.append(fp)

        case _:
            raise ValueError("datasource can be SRTM_GL3 or NASADEM")

    return list(set(files))


def get_dem(
    cells: gpd.GeoDataFrame,
    land: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame | None:
    """Function to get the Digital Elevation Model (DEM) for a given set of h3 cells.
       If no intersection with land then the process returns None

    Args:
        bbox (tuple): The bounding box of the area
        land (shapely.Polygon): The polygon describing world land. Useful to exclude areas on sea

    Returns:
        dem (gpd.GeoDataFrame): GeoDataFrame with the elevation of points in the area
    """

    from rasterio.merge import merge

    # SRTMGL3
    # NASADEM

    minx, miny, maxx, maxy = cells.total_bounds
    # gs = [g["geometry"] for g in json.loads(cells.dissolve().to_json())["features"]]

    if land is not None:
        poly = shapely.Polygon([[minx, miny], [minx, maxy], [maxx, maxy], [maxx, miny]])
        if land.intersects(poly).sum() <= 0:
            return None

    minx = int(np.floor(minx))
    miny = int(np.floor(miny))
    maxx = int(np.ceil(maxx))
    maxy = int(np.ceil(maxy))

    if minx == maxx:
        maxx += 1

    if miny == maxy:
        maxy += 1

    bbox = (minx, miny, maxx, maxy)

    files = download_dem_data(bbox)
    mosaic = []
    meta = {}
    for f in files:
        if not f.exists():
            continue

        if not meta:
            with rasterio.open(f) as raster:
                meta = raster.meta.copy()

        mosaic.append(f)

    if mosaic:
        raster, transform = merge(mosaic, bounds=tuple(cells.total_bounds))
    else:
        return None

    meta.update(
        {
            "driver": "GTiff",
            "height": raster.shape[1],
            "width": raster.shape[2],
            "transform": transform,
        }
    )

    dem = gpd.GeoDataFrame()
    with tempfile.NamedTemporaryFile(suffix=".tif") as f:
        with rasterio.open(f.name, "w", **meta) as r:
            r.write(raster)

        # useful to reduce processing speed (lowers the resolution of the dem)
        resample_factor = 1 / 1

        if resample_factor < 1:
            with rasterio.open(f.name, **meta) as r:
                # resample data to target shape
                raster = r.read(
                    out_shape=(
                        r.count,
                        int(r.height * resample_factor),
                        int(r.width * resample_factor),
                    ),
                    resampling=Resampling.bilinear,
                )

                # scale image transform
                transform = r.transform * r.transform.scale(
                    (r.width / raster.shape[-1]), (r.height / raster.shape[-2])
                )

            # reduce the dem resolution
            meta.update(
                {
                    "driver": "GTiff",
                    "height": raster.shape[1],
                    "width": raster.shape[2],
                    "transform": transform,
                }
            )

            with rasterio.open(f.name, "w", **meta) as r:
                r.write(raster)

        with rioxarray.open_rasterio(f.name, parse_coordinates=True) as r:  # type: ignore
            # clip tot he extents of the cells
            # r = r.rio.clip(gs)
            dem = r.sel(band=1).to_dataframe(name="z")  # type: ignore
            dem["raster_shape"] = [r.sel(band=1).shape] * dem.shape[0]

            dem = dem.reset_index().drop_duplicates(subset=["x", "y"])
            geoms = gpd.points_from_xy(dem["x"], dem["y"])
            dem = gpd.GeoDataFrame(dem, geometry=geoms, crs="epsg:4326").drop(columns=["x", "y"])  # type: ignore
            dem.index.name = "poid"

            dem = gpd.sjoin(
                dem, cells[["geometry"]], predicate="within", how="inner"
            ).drop(columns=["index_right"])

            # TODO must find a better way to pass the info
            dem["raster_transform"] = [transform] * dem.shape[0]

        # calculate the aspect
        with tempfile.NamedTemporaryFile(suffix="-aspect.tif") as ff:
            gdal.DEMProcessing(ff.name, f.name, "aspect")
            with rasterio.open(ff.name, **meta) as r:
                aspect = r.read(1)
                aspect = np.where(aspect == r.nodata, 0, aspect)

                labels = ["n", "e", "s", "w"]
                ln = len(labels)
                labels += [labels[0]]
                bins = [0] + [(180 / ln) + i * (360 / ln) for i in range(ln)] + [360]
                aspect = np.digitize(aspect, bins)
                # last and first group is the same (north)
                aspect = np.where(aspect == aspect.max(), aspect.min(), aspect)
                # fix the dimensions
                dt = rasterio.int16
                aspect = aspect.astype(dt)

                # aspect = rasterio.features.sieve(aspect, size=30, connectivity=4)

                shapes = rasterio.features.shapes(
                    aspect,
                    mask=r.dataset_mask(),
                    transform=r.transform,
                    connectivity=4,
                )
                zones = []
                for s, v in shapes:
                    s = {"type": "Feature", "geometry": s, "properties": {"aspect": v}}
                    zones.append(s)

                if zones:
                    zones = gpd.GeoDataFrame.from_features(zones, crs="epsg:4326")
                else:
                    return None

                zones["aspect"] = (
                    zones["aspect"]
                    .replace({1: "n", 2: "e", 3: "s", 4: "w"})  # type:ignore
                    .astype("category")  # type:ignore
                )

    if dem.empty:
        return None

    dem = (
        gpd.sjoin(
            dem,
            zones,
            how="inner",
            predicate="within",
        )
        .drop(columns=["band", "spatial_ref", "index_right"])
        .drop_duplicates(subset=["geometry"])
    )

    return gpd.GeoDataFrame(dem)


def get_land(fp: str) -> gpd.GeoDataFrame:
    """Function to bring the land boundary of Earth

    Args:
        fp (str): Path to the file containing the world land boundaries

    Returns:
        gpd.Geopandas: A GeoGeoDataFrame with the world's land polygons
    """
    try:
        land = gpd.read_file(fp, engine="pyogrio")
    except:
        import requests
        from fiona.io import ZipMemoryFile

        url = "https://osmdata.openstreetmap.de/download/simplified-land-polygons-complete-3857.zip"
        r = requests.get(url)

        with ZipMemoryFile(r.content) as zip:
            with zip.open(
                "simplified-land-polygons-complete-3857/simplified_land_polygons.shp"
            ) as collection:
                crs = collection.crs
                land = gpd.GeoDataFrame.from_features(collection, crs=crs).to_crs(
                    epsg=4326
                )

                fldr_out = Path(fp).parent
                fldr_out.mkdir(parents=True, exist_ok=True)

                land.to_file(fp, engine="pyogrio")  # type:ignore

    if land.crs.to_epsg() != 4326:  # type:ignore
        land = land.to_crs(epsg=4326)  # type:ignore
    return land  # type:ignore


def get_coastline(fp: str) -> gpd.GeoDataFrame:
    """Function to return the coastline of Earth as described in OSM

    Args:
        fp (str): Path to the file containing the world land boundaries

    Returns:
        gpd.Geopandas: A gGeoDataFrame with the world's coastline
    """
    try:
        coastline = gpd.read_file(fp, engine="pyogrio")
    except:
        import requests
        from fiona.io import ZipMemoryFile

        url = "https://osmdata.openstreetmap.de/download/coastlines-split-4326.zip"
        r = requests.get(url)

        with ZipMemoryFile(r.content) as zip:
            with zip.open("coastlines-split-4326/lines.shp") as collection:
                crs = collection.crs
                coastline = gpd.GeoDataFrame.from_features(collection, crs=crs).to_crs(
                    epsg=4326
                )

                fldr_out = Path(fp).parent
                fldr_out.mkdir(parents=True, exist_ok=True)

                coastline.to_file(fp, engine="pyogrio")  # type:ignore

    if coastline.crs.to_epsg() != 4326:  # type:ignore
        coastline = coastline.to_crs(epsg=4326)  # type:ignore
    return coastline  # type:ignore


def combine_export(
    in_folder: str | Path, out_folder: str | Path, duration: timedelta | None = None
) -> Tuple[int, int, int] | None:
    in_folder = Path(in_folder)

    out_folder = Path(out_folder)
    out_folder.mkdir(parents=True, exist_ok=True)

    ncells = 0
    nzones = 0
    nstations = 0

    cellsf = out_folder / "cells.gpkg"
    if cellsf.exists():
        cellsf.unlink()

    zonesf = out_folder / "zones.gpkg"
    if zonesf.exists():
        zonesf.unlink()

    stationsf = out_folder / "stations.gpkg"
    if stationsf.exists():
        stationsf.unlink()

    for f in list(in_folder.rglob("cells*.gpkg")):
        cells = gpd.read_file(f, engine="pyogrio")
        ncells += cells.shape[0]
        if not cellsf.exists():
            cells.to_file(cellsf, engine="pyogrio")  # type:ignore
        else:
            cells.to_file(cellsf, engine="pyogrio", mode="a")  # type:ignore

    for f in list(in_folder.rglob("zones*.gpkg")):
        zones = gpd.read_file(f, engine="pyogrio")
        nzones += zones.shape[0]
        if not zonesf.exists():
            zones.to_file(zonesf, engine="pyogrio")  # type:ignore
        else:
            zones.to_file(zonesf, engine="pyogrio", mode="a")  # type:ignore

    zones = gpd.read_file(zonesf, engine="pyogrio")
    zones["zone_id"] = zones.reset_index(drop=True).index + 1

    for f in list(in_folder.rglob("stations*.gpkg")):
        stations = gpd.read_file(f, engine="pyogrio")
        nstations += stations.shape[0]
        if not stationsf.exists():
            stations.to_file(stationsf, engine="pyogrio")  # type:ignore
        else:
            stations.to_file(stationsf, engine="pyogrio", mode="a")  # type:ignore

    stations = gpd.read_file(stationsf, engine="pyogrio")
    stations["station_id"] = stations.reset_index(drop=True).index + 1

    if duration:
        _duration = duration.total_seconds()
    else:
        _duration = -1

    pd.Series(
        [ncells, nzones, nstations, _duration],
        index=["cells", "zones", "stations", "duration_s"],
    ).to_csv(out_folder / "summary.csv", header=False)

    return ncells, nzones, nstations
