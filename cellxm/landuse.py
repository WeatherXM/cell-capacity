import geopandas as gpd
import pandas as pd
import shapely

import typing


def prepare_osm(
    osm, min_area_urban: float, min_area_green: float, include_buildings: bool
) -> gpd.GeoDataFrame:
    # quick initial filter based on area. The majority of osm elems are water (not very useful)
    # We exclude those water geoemtries with an area less than 20 times the mimimum green area
    # Lambert azimuthal equal-area = epsg 9822. https://epsg.io/9822-method

    if not include_buildings:
        osm = osm[osm["building"].isna()]  # type:ignore

    osm.geometry = osm.geometry.make_valid()

    osm["type"] = "urban"
    m_water = osm["natural"].isin(["water"])
    osm.loc[m_water, "type"] = "water"
    m_green = osm["leisure"].isin(["park"])
    osm.loc[m_green, "type"] = "green"

    osm = osm.to_crs(epsg=9822)

    m = osm["type"] == "water"  # type:ignore
    water = osm[m]  # type:ignore
    osm = osm[~m]  # type:ignore

    water = water.explode(index_parts=False)  # type:ignore
    water = water[  # type:ignore
        water.area  # type:ignore
        >= 20 * min_area_green
    ]
    water["type"] = "water"

    m = osm["type"] == "urban"  # type:ignore
    urban = osm[m]  # type:ignore
    osm = osm[~m]  # type:ignore

    urban.geometry = urban.buffer(30)
    urban = urban.dissolve().explode(index_parts=False)  # type:ignore
    urban = urban[  # type:ignore
        urban.area  # type:ignore
        >= min_area_urban
    ]
    urban["type"] = "urban"

    m = osm["type"] == "green"  # type:ignore
    green = osm[m]  # type:ignore
    osm = osm[~m]  # type:ignore

    green.geometry = green.buffer(30)
    green = green.dissolve().explode(index_parts=False)  # type:ignore
    green = green[  # type:ignore
        green.area  # type:ignore
        >= min_area_green
    ]
    green["type"] = "green"

    osm = gpd.GeoDataFrame(
        pd.concat([urban, green, water]).reset_index(drop=True)  # type:ignore
    ).to_crs(epsg=4326)

    cols = ["geometry", "type"]
    osm = osm[cols]  # type:ignore

    osm.geometry = osm.geometry.make_valid()
    mask = osm.geom_type.isin(["Polygon", "MultiPolygon"])  # type:ignore
    osm = osm.loc[mask]  # type:ignore
    osm.geometry = shapely.set_precision(osm.geometry.values, 1e-9)  # type:ignore

    return osm


@typing.no_type_check
def landuse_from_osm(
    osm: gpd.GeoDataFrame,
    land: gpd.GeoDataFrame,
    coastline: gpd.GeoDataFrame,
    cells: gpd.GeoDataFrame,
    min_area_urban: float,
    min_area_green: float,
    largest_urban: bool,
    largest_green: bool,
) -> gpd.GeoDataFrame | None:
    """Create landuse polygons from OSM (Openstreetmap) data.
       The idea is to use all urban-like, green-like (e.g. parks) and coastline elements form OSM and classify their landuse as urban, green and coastline respectively.
       The remaining area of the cell is classified as landuse: rural

       The query which returned the OSM data is:
       landuse IN ('residential', 'commercial', 'retail') OR
       place IN ('neighbourhood') OR aeroway IN ('aerodrome') OR
       leisure IN ('park') OR natural IN ('water')

       Although not fully inclusive, it has proven to achieve adequate
       coverage of the landuse elements we are interested in

    Args:
        osm_elems (gpd.GeoDataFrame): OSM elements (created in previous step)
        land (gpd.GeoDataFrame): The land boundary
        coastline (gpd.GeoDataFrame): The coastline boundary
        cells (gpd.GeoDataFrame): H3 cells to run the analysis
        min_area_urban (float): The minimum area to consider an urban area as individual element
        min_area_green (float): The minimum area to consider a green area as individual element
        largest_urban (bool): Boolean to specify which urban areas to keep. If True only the largest urban area for each cell will be exported else all
        largest_green (bool): Boolean to specify which green areas to keep. If True only the largest green area for each urban region will be exported else all

    Returns:
        gpd.GeoDataFrame: GeoDataFrame with all the different landuses in the analysed area
    """

    def _keep_largest_landuse(landuse, landuse_type):
        # keep only the largest landuse area per type and per cell

        if "h3_id" in landuse:
            del landuse["h3_id"]

        m = landuse["type"] == landuse_type
        lu_type = landuse[m]
        landuse = landuse[~m]

        lu_type_copy = lu_type.copy()
        lu_type_copy.geometry = lu_type_copy.to_crs(crs).geometry.centroid.to_crs(
            epsg=4326
        )
        lu_type_copy = lu_type_copy.sjoin(
            cells[["geometry"]], predicate="intersects"
        ).rename(columns={"index_right": "h3_id"})
        lu_type = gpd.GeoDataFrame(lu_type.join(lu_type_copy[["h3_id"]]))
        idxs = lu_type.groupby("h3_id", observed=True)["area"].idxmax().values
        lu_type = lu_type.loc[idxs]  # type:ignore
        landuse = pd.concat([landuse, lu_type]).reset_index(drop=True)
        return landuse

    crs = cells.estimate_utm_crs()

    urban_elems = osm.loc[osm["type"] == "urban"]  # type:ignore
    if urban_elems.empty:
        landuse = cells.dissolve("h3_id", observed=True).reset_index()
        landuse["type"] = "rural"
        landuse["area"] = landuse.to_crs(crs).area  # type:ignore
        landuse["eligible"] = True
        landuse.index.name = "luid"
        return gpd.GeoDataFrame(landuse)

    green_elems = osm.loc[osm["type"] == "green"]  # type:ignore

    urban_elems = urban_elems.to_crs(crs)  # type: ignore
    urban_elems.loc[:, "landuse_type"] = "urban"  # type:ignore

    water_elems = osm.loc[osm["type"] == "water"]  # type:ignore
    if not water_elems.empty:
        water_elems = gpd.GeoDataFrame(water_elems)
        water_elems["landuse_type"] = "water"
        try:
            land = land.overlay(
                water_elems[["geometry"]].to_crs(epsg=4326),  # type:ignore
                how="difference",
                keep_geom_type=True,
            )  # type:ignore
        except:
            pass

    # join nearby areas
    urban_elems.geometry = urban_elems.buffer(30)  # type:ignore
    urban_elems = urban_elems.to_crs(epsg=4326)  # type:ignore
    urban_elems = urban_elems.dissolve().explode(ignore_index=True)
    # First pass filtering, to improve performance
    urban_elems = urban_elems[urban_elems.to_crs(crs).area > min_area_urban]
    urban_elems.geometry = urban_elems.geometry.make_valid()

    df = cells.reset_index().overlay(
        land, how="intersection", keep_geom_type=True
    )  # type:ignore
    df.geometry = df.geometry.make_valid()
    df.geometry = shapely.set_precision(df.geometry.values, 1e-9)

    if df.empty:
        return None

    cols = ["h3_id", "geometry"]

    urban = df.overlay(urban_elems, how="intersection", keep_geom_type=True).explode(
        ignore_index=True
    )[cols]

    if not urban.empty:
        urban.loc[:, "landuse_type"] = "urban"

    green = gpd.GeoDataFrame()

    if not urban.empty:
        urban_cells = list(urban["h3_id"].unique())
        green = green_elems.to_crs(crs)  # type:ignore
        green.geometry = green.buffer(20)  # type: ignore
        green = green.dissolve().explode(ignore_index=True).to_crs(epsg=4326)  # type: ignore
        green = df.overlay(green, how="intersection", keep_geom_type=True)[cols]
        green = green.reset_index(drop=True)

        if not green.empty:
            green.loc[:, "landuse_type"] = "green"

        green = green[green["h3_id"].isin(urban_cells)]

    # rural is the difference between the cell and the union of urban and the water
    # other = pd.concat([urban, water])
    rural = df.overlay(urban_elems, how="difference", keep_geom_type=True)[cols]
    rural.loc[:, "landuse_type"] = "rural"

    # include coastline only for urban areas
    if not urban.empty:
        coastline.geometry = coastline.to_crs(crs).buffer(50)  # type:ignore
        coastline = cells.reset_index().overlay(  # type:ignore
            coastline.to_crs(epsg=4326), how="intersection", keep_geom_type=True
        )

        if not coastline.empty:
            coastline = coastline.rename(
                columns={"index_right": "h3_id"}
            )  # type:ignore
            coastline.loc[:, "landuse_type"] = "coastline"
        else:
            coastline = gpd.GeoDataFrame()
    else:
        coastline = gpd.GeoDataFrame()

    # explode and calculate area
    landuse = gpd.GeoDataFrame(
        pd.concat([urban, rural, green, coastline]).rename(
            columns={"landuse_type": "type"}
        )  # type:ignore
    ).explode(ignore_index=True)

    # I encountered error with invalid geometry
    landuse.geometry = shapely.set_precision(landuse.geometry.values, 1e-9)
    landuse.geometry = landuse.make_valid()  # type: ignore
    landuse = (
        landuse.dissolve(["type"])  # type:ignore
        .dropna(subset=["geometry"])
        .explode(index_parts=False)
        .reset_index()
        .drop(columns=["h3_id"])  # type:ignore
    )  # type:ignore

    mask = landuse.geom_type.isin(["Polygon", "MultiPolygon"])
    landuse = landuse.loc[mask]  # type:ignore

    # cut landuse by the cells
    landuse = cells.overlay(  # type:ignore
        landuse,
        how="intersection",
        keep_geom_type=True,
    ).explode(
        ignore_index=True  # type:ignore
    )

    landuse["area"] = landuse.to_crs(
        crs
    ).area  # type:ignore # do not drop, required later

    if min_area_urban:
        mask = (landuse["type"] == "urban") & (landuse["area"] < min_area_urban)
        landuse.loc[mask, "type"] = None  # type:ignore
        landuse = landuse.dropna(subset=["type"])

    if min_area_green:
        mask = (landuse["type"] == "green") & (landuse["area"] < min_area_green)
        landuse.loc[mask, "type"] = None  # type:ignore
        landuse = landuse.dropna(subset=["type"])

    # Remove green and coastline if no urban area in each cell left after filtering
    sjoin = gpd.sjoin(
        landuse.representative_point().to_frame(name="geometry"),  # type:ignore
        cells,
        predicate="within",
    )
    for _, g in sjoin.groupby("index_right"):
        idxs = g.index
        lu = landuse.loc[idxs]
        if lu[lu["type"] == "urban"].empty:
            # drop the coastlines
            _cs = lu[lu["type"] == "coastline"].index
            landuse = landuse.drop(index=_cs)

            _gs = lu[lu["type"] == "green"].index
            # convert green to rural
            landuse.loc[_gs, "type"] = "rural"

    if largest_urban:
        landuse = _keep_largest_landuse(landuse, "urban")

    largest_coastline = True
    if largest_coastline:
        landuse = _keep_largest_landuse(landuse, "coastline")

    # Largest green differs because is it relative to each
    if largest_green:
        landuse = _keep_largest_landuse(landuse, "green")

    c = "index_right"
    if c in landuse:
        landuse = landuse.drop(columns=[c])

    landuse["eligible"] = True
    mask = landuse["type"].isin(["water"])
    landuse.loc[mask, "eligible"] = False
    landuse.index.name = "luid"

    return landuse  # type:ignore
