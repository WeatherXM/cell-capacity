import json
from typing import Tuple
import geopandas as gpd
import h3
import pandas as pd
import shapely


def assign_h3_cell_to_dem(
    dem: gpd.GeoDataFrame, cells: gpd.GeoDataFrame
) -> gpd.GeoDataFrame:
    """Assign the cell id for each point in the Digital Elevation Model (dem) GeoDataFrame
       If cell_buffer is greater than 0, then the dem points inside the buffer are filtered out
    Args:
        dem (gpd.GeoDataFrame): The digital elevation model
        cells (gpd.GeoDataFrame): The h3 cells
        cell_buffer (float): The buffer area inside the cell to avoid placing stations
    Returns:
        gpd.GeoDataFrame: DEM with collumn h3_id
    """
    dem = (
        gpd.sjoin(dem.reset_index(), cells, predicate="within", how="inner")
        .rename(columns={"index_right": "h3_id"})
        .reset_index(drop=True)
    ).set_index("poid")

    return dem


def h3_cells_from_bbox(
    bbox: Tuple[int, int, int, int],
    resolution: int,
) -> gpd.GeoDataFrame:
    """Create the GeoDataFrame of the h3 cells in the given bounding box
    Args:
        bbox (Tuple[int, int, int, int]): Tuple of the bounding box in the format (min longitude , min latitude , max longitude , max latitude)
        resolution (int, optional): The resilution of the h3 grid. Defaults to 7.
    Returns:
        gpd.GeoDataFrame: The GeoDataFrame containing the h3 cells
    """
    minx, miny, maxx, maxy = bbox
    poly = h3.Polygon([[miny, minx], [miny, maxx], [maxy, maxx], [maxy, minx]])

    cells = []
    for cell_id in h3.polygon_to_cells(poly, res=resolution):
        c = shapely.Polygon(h3.cell_to_boundary(cell_id, geo_json=True))
        cells.append((cell_id, c))

    cells = pd.DataFrame(cells, columns=["h3_id", "geom"])
    cells = gpd.GeoDataFrame(
        cells, geometry=cells["geom"], crs="epsg:4326"
    ).drop(  # type: ignore
        columns=["geom"]
    )

    cells["h3_id"] = cells["h3_id"].astype("category")  # type: ignore
    cells = cells.set_index("h3_id")
    cells.geometry = cells.geometry.make_valid()
    cells.geometry = shapely.set_precision(cells.geometry.values, 1e-9)

    return cells  # type: ignore
