from dotenv import load_dotenv

load_dotenv()

import os
import sys

if not sys.warnoptions:
    import warnings

    warnings.simplefilter("ignore")
    os.environ["PYTHONWARNINGS"] = "ignore"


import datetime
import pandas as pd
import geopandas as gpd
import itertools
import logging
import more_itertools
import numpy as np
import pathlib
import tempfile
import toml
import typer
import uuid

# from multiprocessing import Pool
from multiprocessing import get_context

from rich.console import Console
from rich.table import Table
from rich.progress import (
    Progress,
    SpinnerColumn,
    MofNCompleteColumn,
    TimeElapsedColumn,
    TextColumn,
)

from typing import Tuple, Optional
from typing_extensions import Annotated

from cellxm.io import (
    download_osm_elems,
    get_dem,
    get_land,
    get_coastline,
    combine_export,
)
from cellxm.locators import (
    locate_station_zones,
)
from cellxm.h3_helpers import assign_h3_cell_to_dem, h3_cells_from_bbox
from cellxm.landuse import prepare_osm, landuse_from_osm
from cellxm.utils import categoricals_to_str, overwrite_config

debug = os.getenv("DEBUG", "False").lower() in ("true", "1", "t")

app = typer.Typer()


def locate_stations_wrap(
    cells: gpd.GeoDataFrame,
    land: gpd.GeoDataFrame,
    coastline: gpd.GeoDataFrame,
    osm: gpd.GeoDataFrame,
    elev_diff_thresh: float,
    cell_buffer: float,
    min_area_green: float,
    min_area_urban: float,
    largest_urban: bool,
    largest_green: bool,
    min_aspect_zone_area_cell_perc: float,
) -> Optional[Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame]]:
    """Wrapper function to assign to multiprocessing"""

    try:
        cellsbuf = cells.copy()
        if cell_buffer:
            crs = cells.estimate_utm_crs()
            cellsbuf.geometry = (
                cells.to_crs(crs).buffer(-cell_buffer).to_crs(epsg=4326)  # type:ignore
            )

        cellsbuf = gpd.GeoDataFrame(cellsbuf)

        dem = get_dem(cells=cellsbuf, land=land)

        if dem is None:
            return None

        dem = assign_h3_cell_to_dem(dem, cells)

        if dem.empty:
            return None

        landuse = landuse_from_osm(
            osm,
            land,
            coastline,
            cellsbuf,
            min_area_green=min_area_green,
            min_area_urban=min_area_urban,
            largest_urban=largest_urban,
            largest_green=largest_green,
        )

        if landuse is None:
            return None

        zones = locate_station_zones(
            dem, cells, landuse, elev_diff_thresh, min_aspect_zone_area_cell_perc
        )

        return cells, zones
    except BaseException as e:
        # print(e)
        print(f"error at bbox: {cells.total_bounds}")
        return None


def call_locate_stations_wrap(params):
    return locate_stations_wrap(*params)


if not debug:

    @app.command()
    def combine(in_folder: pathlib.Path, out_folder: pathlib.Path):
        r = combine_export(in_folder, out_folder, duration=None)
        if r:
            ncells, nzones, nstations = r
            print(
                f"combined and exported input of folder {in_folder} to folder: {out_folder}"
            )
            print(f"n-cells: {ncells}, n-zones: {nzones}, n-stations: {nstations}")
        else:
            logging.warn(f"No results found in {in_folder}")


@app.command()
def locate(
    country_code: str = typer.Argument(
        "GR",
        help="The country code to run the application. Should follow the ISO 3166-1 alpha-2 specification",
    ),
    config: str = typer.Argument(
        "./config.toml", help="The location of the config file"
    ),
    h3_resolution: int = typer.Option(
        None, help="The resolution for the H3 grid. Should be between 0 and 15"
    ),
    elev_diff_thresh: float = typer.Option(
        None, help="The elevation difference inside a cell to trigger the heuristic"
    ),
    cell_buffer: float = typer.Option(
        None,
        help="Buffer distance around a cell where stations cannot be placed.",
    ),
    min_area_green: float = typer.Option(
        None, help="The minimum area to consider a green area as individual element"
    ),
    min_area_urban: float = typer.Option(
        None, help="The minimum area to consider an urban area as individual element"
    ),
    min_aspect_zone_area_cell_perc: float = typer.Option(
        None,
        help="The minimum area of an aspect zone as a percentage of the average h3 cell area",
    ),
    largest_urban: Annotated[
        bool,
        typer.Option(
            "--largest-green",
            help="Boolean to specify which urban areas to keep. If True only the largest urban area for each cell will be exported else all",
        ),
    ] = False,
    largest_green: Annotated[
        bool,
        typer.Option(
            "--largest-urban",
            help="Boolean to specify which green areas to keep. If True only the largest green area for each urban region will be exported else all",
        ),
    ] = False,
    include_buildings: bool = typer.Option(
        None,
        help="""Use buildings to infer urban areas. This reduces the likelihood of not accurately capturing urban areas due to missing information in OSM.
                       However this comes at increased computational cost (~1.5 longer processing time) """,
    ),
    out_folder: str = typer.Option(
        "./data/outputs", help="Location to export the results"
    ),
    ncpus: int = typer.Option(None, help="The number of cpus to be used"),
) -> Optional[Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]]:
    """The main function to locate weather stations on a country. It follows the methodology described in the Methodology Report

    Args:
        country_code (str): The country code to run the application. Should follow the ISO 3166-1 alpha-2 specification. Defaults to "GR".
        config (str, optional): The location of the config file. Defaults to "./config.toml".
        h3_resolution (int): The resolution for the H3 grid. Should be between 0 and 15.
        elev_diff_thresh (float): The elevation difference inside a cell to trigger the heuristic.
        cell_buffer (float): Buffer distance around a cell where stations cannot be placed.
        min_area_green (float): The minimum area to consider a green area as individual element.
        min_area_urban (float): The minimum area to consider an urban area as individual element.
        min_aspect_zone_area_cell_perc (float): The minimum area of an aspect zone as a percentage of the average h3 cell area.
        largest_urban (bool): Boolean to specify which urban areas to keep. If True only the largest urban area for each cell will be exported else all. Defaults to False
        largest_green (bool): Boolean to specify which green areas to keep. If True only the largest green area for each urban region will be exported else all. Defaults to False
        out_folder (str): Location to export the results. Defaults to "./data/outputs".
        ncpus (int): The number of cpus to be used. Defaults to None

    Returns:
        Optional[Tuple[gpd.GeoDataFrame, gpd.GeoDataFrame, gpd.GeoDataFrame]]: A tuple of cells, zones and stations
    """
    with open(config, "r") as f:
        configd = toml.load(f)

    configd = overwrite_config(
        configd,
        h3_resolution,
        elev_diff_thresh,
        cell_buffer,
        min_area_urban,
        min_area_green,
        largest_urban,
        largest_green,
        min_aspect_zone_area_cell_perc,
        include_buildings,
        out_folder,
    )

    h3_resolution = configd["h3_resolution"]
    elev_diff_thresh = configd["elev_diff_thresh"]
    cell_buffer = configd["cell_buffer"]
    min_area_urban = configd["min_area_urban"]
    min_area_green = configd["min_area_green"]
    largest_urban = configd["largest_urban"]
    largest_green = configd["largest_green"]
    min_aspect_zone_area_cell_perc = configd["min_aspect_zone_area_cell_perc"]
    include_buildings = configd["include_buildings"]
    out_folder = configd["outfolder"]
    ncpus = configd.get("ncpus", None)

    if not ncpus:
        ncpus = os.cpu_count() or 1

    configd["ncpus"] = ncpus

    logging.info(f"executing script with params:{configd.values}")

    # if debug:
    #     # for testing
    #     country_code = "AD"
    #     h3_resolution = 7
    #     elev_diff_thresh = 100
    #     # snap_distance = 200  # type:ignore
    #     min_area_green = 15000
    #     min_area_urban = 400000
    #     cell_buffer = 250
    #     largest_urban = True
    #     largest_green = True
    #     # min_aspect_zone_area_cell_perc = 0.16
    #     out_folder = f"./data/outputs/{country_code}-test"
    #     # ncpus = 1

    if not (0 <= h3_resolution <= 15):
        raise ValueError("h3_resolution should be between 0 and 15")

    if not (0 < min_aspect_zone_area_cell_perc <= 1):
        raise ValueError("min_aspect_zone_area_cell_perc must be between 0 and 1")

    country_code = country_code.upper()

    cachedir = pathlib.Path("./.cache")
    procdir = cachedir / "processing"
    procdir.mkdir(parents=True, exist_ok=True)

    start = datetime.datetime.now()
    with Progress(
        SpinnerColumn(),
        TextColumn("[progress.description]{task.description}"),
        TimeElapsedColumn(),
    ) as progress:
        task_id = progress.add_task("[cyan]Preprocessing...")

        boundary = download_osm_elems(country_code)

        minx, miny, maxx, maxy = boundary.bounds
        land = get_land(fp="./data/inputs/land.gpkg")
        coast = get_coastline(fp="./data/inputs/coastline.gpkg")

        land = land.clip(boundary, keep_geom_type=True)
        coast = coast.clip(boundary, keep_geom_type=True)

        step = 1
        rngx = list(range(int(np.floor(minx)), int(np.ceil(maxx))))  # type:ignore
        rngy = list(range(int(np.floor(miny)), int(np.ceil(maxy))))  # type:ignore

    with tempfile.TemporaryDirectory(dir=procdir) as tdir:
        tdir = pathlib.Path(tdir)
        gxys = list(itertools.product(rngx, rngy))

        osms = []
        for x, y in gxys:
            # print(x, y)
            bbox = (x, y, x + step, y + step)
            osmc = gpd.read_file(cachedir / f"osm/{country_code}.gpkg", bbox=bbox)
            osmc = prepare_osm(osmc, min_area_urban, min_area_green, include_buildings)
            osms.append(osmc)

        if osms:
            osm = pd.concat(osms).reset_index(drop=True)
        else:
            osm = gpd.GeoDataFrame()

        with Progress(
            SpinnerColumn(),
            *Progress.get_default_columns(),
            MofNCompleteColumn(),
            TimeElapsedColumn(),
        ) as progress:
            task_id = progress.add_task(
                "[cyan]Identifying stations...", total=len(gxys)
            )

            for xys in more_itertools.batched(gxys, ncpus * 2):
                args = []
                for xy in xys:
                    x, y = xy

                    # At edge cases (e.g. around the poles) cells can cause issues
                    if (abs(x) >= 178) or (abs(y) >= 88):
                        progress.update(task_id, advance=1)
                        continue

                    bbox = (x, y, x + step, y + step)

                    cells = h3_cells_from_bbox(bbox, resolution=h3_resolution)
                    cells = gpd.GeoDataFrame(cells[cells.centroid.within(boundary)])  # type: ignore

                    if cells.empty:
                        progress.update(task_id, advance=1)
                        continue

                    landc = land.clip(cells, keep_geom_type=True)  # type:ignore
                    if landc.empty:
                        progress.update(task_id, advance=1)
                        continue

                    coastc = coast.clip(cells, keep_geom_type=True)  # type:ignore
                    osmc = osm.clip(cells, keep_geom_type=True)  # type:ignore

                    args.append(
                        [
                            cells,
                            landc,
                            coastc,
                            osmc,
                            elev_diff_thresh,
                            cell_buffer,
                            min_area_green,
                            min_area_urban,
                            largest_urban,
                            largest_green,
                            min_aspect_zone_area_cell_perc,
                        ]
                    )

                if not args:
                    continue

                # fix weird freezing
                with get_context("spawn").Pool(ncpus) as pool:
                    for result in pool.imap_unordered(
                        call_locate_stations_wrap,
                        args,
                    ):
                        if result:
                            cells, zones = result

                            if cells.empty:
                                continue

                            if zones.empty:
                                continue

                            rnd = uuid.uuid4()

                            cells.geometry = cells.geometry.make_valid()
                            cells["country"] = country_code

                            # keep zones whose centroid is within the boundary
                            zones.geometry = zones.geometry.make_valid()
                            m = zones.representative_point().within(boundary)
                            zones = zones[m]

                            m = zones.geom_type.isin(  # type:ignore
                                ["Polygon", "MultiPolygon"]
                            )
                            zones = zones[m]  # & m2]  # type:ignore
                            zones = zones.reset_index(drop=True)  # type:ignore
                            zones.index.name = "zone_id"  # type:ignore
                            zones["country"] = country_code
                            zones = gpd.GeoDataFrame(zones.reset_index())

                            h3_ids = list(zones["h3_id"].unique())  # type:ignore
                            if h3_ids:
                                cells = cells.loc[h3_ids]

                            stations = zones.copy()
                            stations.geometry = (
                                stations.representative_point()
                            )  # type:ignore
                            stations.index.name = "station_id"
                            stations = gpd.GeoDataFrame(stations)

                            cells = categoricals_to_str(
                                cells.reset_index()  # type:ignore
                            )  # type:ignore
                            zones = categoricals_to_str(zones)
                            stations = categoricals_to_str(stations)

                            cells.to_file(  # type:ignore
                                tdir / f"cells-{rnd}.gpkg",
                                engine="pyogrio",
                                index=False,
                            )

                            zones.to_file(  # type:ignore
                                tdir / f"zones-{rnd}.gpkg",
                                engine="pyogrio",
                                index=False,  # type:ignore
                            )

                            stations.to_file(  # type:ignore
                                tdir / f"stations-{rnd}.gpkg",
                                engine="pyogrio",
                                index=False,  # type:ignore
                            )

                        progress.update(task_id, advance=1)

        progress.add_task("[cyan]Exporting...")

        end = datetime.datetime.now()
        duration = end - start
        r = combine_export(tdir, out_folder, duration)

        if r:
            with open(pathlib.Path(out_folder) / "config.toml", "w") as f:
                toml.dump(configd, f)

            ncells, nzones, nstations = r
        else:
            logging.info("no stations were identified")
            return

    table = Table(title="Results")
    table.add_column("Result type", justify="left", style="cyan")
    table.add_column("#", style="green")

    table.add_row("cells", f"{ncells}")
    table.add_row("zones/stations", f"{nstations}")

    console = Console()
    console.print(table)
    return ncells, nzones, nstations  # type: ignore


if __name__ == "__main__":
    app()
