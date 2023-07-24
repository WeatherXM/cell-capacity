from typing import Tuple

import numpy as np
import pandas as pd
import geopandas as gpd

# from cellxm import ogr2ogr

import rasterio.features
from rasterio.enums import MergeAlg

from shapely import Polygon


def close_holes(poly: Polygon) -> Polygon:
    """
    Close polygon holes by limitation to the exterior ring.
    Args:
        poly: Input shapely Polygon
    Example:
        df.geometry.apply(lambda p: close_holes(p))
    """
    if poly.interiors:
        return Polygon(list(poly.exterior.coords))
    else:
        return poly


def locate_station_zones(
    dem: gpd.GeoDataFrame,
    cells: gpd.GeoDataFrame,
    landuse: gpd.GeoDataFrame,
    elev_diff_thresh: float,
    min_aspect_zone_area_cell_perc: float,
) -> gpd.GeoDataFrame:
    """Function to locate the zones in which station will be placed

    Args:
        dem (gpd.GeoDataFrame): DEM
        landuse (gpd.GeoDataFrame): landuse info
    """

    zones = []

    # keep only the areas where a station can be placed
    lu = landuse[landuse["eligible"]]  # type:ignore
    mask_topo = lu["type"] == "rural"  # topographic analysis # type:ignore
    mask_urban = ~mask_topo  # Urban landuse analysis

    topo = gpd.GeoDataFrame(lu[mask_topo])  # type:ignore
    urban = gpd.GeoDataFrame(lu[mask_urban])  # type:ignore

    crs = cells.estimate_utm_crs()

    avg_cell_area = cells.to_crs(crs).area.mean()  # type:ignore
    min_aspect_zone_area = min_aspect_zone_area_cell_perc * avg_cell_area

    if not topo.empty:
        demf = (
            gpd.sjoin(dem, topo[["geometry"]], how="inner", predicate="within")
            .drop(columns=["index_right"])
            .drop_duplicates()
        )

        if not demf.empty:
            pzones = locate_zones_topography(
                demf, cells, elev_diff_thresh, min_aspect_zone_area
            )

            zones.append(pzones)

    if not urban.empty:
        demf = (
            gpd.sjoin(dem, urban[["geometry"]], how="inner", predicate="within")
            .drop(columns=["index_right"])
            .drop_duplicates()
        )
        szones = locate_zones_urban_landuse(urban)  # type:ignore
        # szones = rasterize_landuse(demf, szones)
        zones.append(szones)

    if zones:
        zones = gpd.GeoDataFrame(pd.concat(zones).reset_index(drop=True)).drop(
            columns=["h3_id"]
        )
    else:
        zones = gpd.GeoDataFrame()

    # # identify_subzones
    m = zones["type"] == "lu-green"
    green = (
        zones[m]  # type:ignore
        .representative_point()  # type:ignore
        .to_frame(name="geometry")
    )
    m = zones["type"].str.startswith("topo-")  # type:ignore
    topo = (
        zones[m]
        .geometry.apply(lambda p: close_holes(p))  # type:ignore
        .to_frame(name="geometry")
    )

    s = gpd.sjoin(green, topo, predicate="within")
    subzone_idxs = s[s.index != s["index_right"]].index

    new_types = zones.loc[s["index_right"].values, "type"].values  # type:ignore
    zones.loc[subzone_idxs, "type"] = new_types  # type:ignore

    zones = rasterize_zones(dem, zones, sieve_size=4)  # type:ignore
    zones = zones.sjoin(cells).rename(columns={"index_right": cells.index.name})

    # rasterization creates some new zones which mst be filtered
    zones["area"] = zones.to_crs(crs).area

    # special case. If a cell contains only topo-h and urban zones do not return the aspect with the lowest points
    # expect if the lowest aspect zone is larger than the threshold
    drps = []
    for _, g in zones.groupby("h3_id", observed=True):
        if g.empty:
            continue

        d = []
        s1 = set(["topo-h", "lu-urban"])
        s2 = set(g["type"].values)
        if s1.issubset(s2):
            d = list(s2.difference(s1))

        if d:
            if "lu-green" in d:
                d.remove("lu-green")

            # drop only when there is one aspect zone left
            if (len(d) != 1) or (not d[0].startswith("topo-")):
                continue

            g = g[g["type"].isin(d)]
            drps.extend(g[g["area"] < min_aspect_zone_area].index.values)

    if drps:
        zones = zones.drop(index=set(drps))

    idxs = zones.groupby(["h3_id", "type"], observed=True)["area"].idxmax()
    zones = zones.loc[idxs].drop(columns=["area"])

    return zones  # type:ignore


def locate_zones_topography(
    dem: gpd.GeoDataFrame,
    cells: gpd.GeoDataFrame,
    elev_diff_thresh: float,
    min_aspect_zone_area: float,
) -> gpd.GeoDataFrame:
    """Locate the station zones based on topographic analysis (see the methodology doc for more info)

    Args:
        dem (gpd.GeoDataFrame): DEM
        cells (gpd.GeoDataFrame): H3 cells
        elev_diff_thresh (float): Minimum elevation difference threshold to trigger topography analysis

    Returns:
        gpd.GeoDataFrame: Station zones from topographic analysis
    """

    crs = cells.estimate_utm_crs()
    cells_highest = dem.groupby("h3_id", observed=True)["z"].max()

    # filter out points with elevation difference less than the threshold
    df = pd.merge(
        dem,
        cells_highest,
        how="left",
        left_on="h3_id",
        right_index=True,
        suffixes=("", "_cell_highest"),
    )

    mask_low = df["z_cell_highest"] - df["z"] >= elev_diff_thresh

    # The high zone covers half the elevation difference
    mask_high = df["z_cell_highest"] - df["z"] < elev_diff_thresh / 2

    rough_cells = list(df.loc[mask_low, "h3_id"].unique())
    # flat_cells = list(set(df["h3_id"]) - set(rough_cells))

    d = {0: "n", 1: "e", 2: "s", 3: "w"}
    d_inv = {v: k for k, v in d.items()}

    # TODO Assumes only one raster size! It also creates problem with split cells at boundaries
    rs = df["raster_shape"].value_counts().idxmax()
    rt = df["raster_transform"].value_counts().idxmax()

    df_low = df[mask_low]
    df_high = df[mask_high]
    df_low = df_low[df_low["h3_id"].isin(rough_cells)]
    # df_high = df_high[df_high["h3_id"].isin(rough_cells)]

    # start with high ground zones
    df_high.loc[:, "aspect"] = 0
    # create tuples of geometry, value pairs, where value is the attribute value you want to burn
    gvs = ((g, v) for g, v in zip(df_high.geometry, df_high["aspect"]))

    # Rasterize vector using the shape and transform of the raster
    rasterized = rasterio.features.rasterize(
        gvs,
        out_shape=rs,
        transform=rt,  # type: ignore
        all_touched=True,
        fill=-1,  # background value
        merge_alg=MergeAlg.replace,
        dtype=rasterio.int16,
        default_value=-1,
    )

    mask = np.where(rasterized == -1, False, True)
    shapes = rasterio.features.shapes(rasterized, mask=mask, transform=rt, connectivity=4)  # type: ignore
    zones = []
    for s, v in shapes:
        s = {"type": "Feature", "geometry": s, "properties": {"aspect": "h"}}
        zones.append(s)

    if not df_low.empty:
        # create tuples of geometry, value pairs, where value is the attribute value you want to burn
        gvs = ((g, v) for g, v in zip(df_low.geometry, df_low["aspect"].replace(d_inv)))

        # Rasterize vector using the shape and transform of the raster
        rasterized = rasterio.features.rasterize(
            gvs,
            out_shape=rs,
            transform=rt,  # type: ignore
            all_touched=True,
            fill=-1,  # background value
            merge_alg=MergeAlg.replace,
            dtype=rasterio.int16,
            default_value=-1,
        )

        mask = np.where(rasterized == -1, False, True)

        shapes = rasterio.features.shapes(rasterized, mask=mask, transform=rt, connectivity=4)  # type: ignore

        for s, v in shapes:
            s = {"type": "Feature", "geometry": s, "properties": {"aspect": d[v]}}
            zones.append(s)

    zones = gpd.GeoDataFrame.from_features(zones).set_crs(dem.crs)
    zones = zones.overlay(
        cells, how="intersection", keep_geom_type=True
    ).explode(  # type:ignore
        ignore_index=True
    )  # type:ignore

    zones["centroid"] = zones.geometry.representative_point()
    zones = zones.set_geometry("centroid")  # type:ignore
    zones = gpd.sjoin(zones, cells, predicate="within").rename(
        columns={"index_right": "h3_id"}
    )
    zones = zones.set_geometry("geometry").drop(columns=["centroid"])
    zones["area"] = zones.to_crs(zones.estimate_utm_crs()).area

    # for each cell keep the aspect with the lowest point
    idxs_lowest = df.groupby("h3_id", observed=True)["z"].idxmin().unique()
    df_lowest = df.loc[idxs_lowest]
    zones_lowest_index = gpd.sjoin(zones, df_lowest, predicate="contains").index

    zones["keep"] = False
    zones.loc[zones_lowest_index, "keep"] = True

    zones = zones[
        (zones["aspect"] == "h")
        | zones["keep"]
        | (zones["area"] >= min_aspect_zone_area)
    ]

    zones = zones.reset_index(drop=True)
    zones.index.name = "szid"

    idxs = zones.groupby(["h3_id", "aspect"], observed=True)["area"].idxmax().unique()
    zones = zones.loc[idxs].drop(columns=["area", "keep"])

    zones.geometry = zones.to_crs(crs).buffer(10).to_crs(epsg=4326)
    zones["type"] = "topo-" + zones["aspect"].astype(str)
    zones = zones.drop(columns=["aspect"])

    return zones  # type: ignore


def locate_zones_urban_landuse(
    landuse: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    """
    Locate the station zones based on urban land use analysis (see the methodology doc for more info)

    Args:
        dem (gpd.GeoDataFrame): _description_
        landuse (gpd.GeoDataFrame): _description_

    Returns:
        gpd.GeoDataFrame: The station zones from urban land use
    """
    zones = landuse.reset_index()  # type:ignore
    cols = ["h3_id", "geometry", "type"]
    zones = zones[cols]
    zones["type"] = "lu-" + zones["type"].astype(str)
    return zones  # type:ignore


def rasterize_zones(
    dem: gpd.GeoDataFrame, zones: gpd.GeoDataFrame, sieve_size: int | None = None
):
    shape = dem["raster_shape"].value_counts().idxmax()  # type: ignore
    transform = dem["raster_transform"].value_counts().idxmax()  # type: ignore

    cat = zones["type"].astype("category").cat  # type: ignore
    d = dict(zip(cat.codes, zones["type"]))  # type: ignore

    gvs = ((g, v) for g, v in zip(zones.geometry, cat.codes))

    # Rasterize vector using the shape and transform of the raster
    rasterized = rasterio.features.rasterize(
        gvs,
        out_shape=shape,
        transform=transform,  # type: ignore
        all_touched=True,
        fill=-1,  # background value
        merge_alg=MergeAlg.replace,
        dtype=rasterio.int16,
        default_value=-1,
    )

    if sieve_size:
        mask = np.where(rasterized == -1, False, True)
        rasterized = rasterio.features.sieve(
            rasterized, size=sieve_size, connectivity=4
        )

    mask = np.where(rasterized == -1, False, True)
    shapes = rasterio.features.shapes(
        rasterized,
        mask=mask,
        transform=transform,  # type: ignore
        connectivity=4,
    )  # type: ignore

    zones = []  # type: ignore
    for s, v in shapes:
        s = {"type": "Feature", "geometry": s, "properties": {"type": d[v]}}
        zones.append(s)  # type: ignore

    zones = gpd.GeoDataFrame.from_features(zones).set_crs(dem.crs)  # type: ignore
    return zones
