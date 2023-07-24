import geopandas as gpd


def overwrite_config(
    config,
    h3_resolution,
    elev_diff_thresh,
    cell_buffer,
    min_area_urban,
    min_area_green,
    largest_urban,
    largest_green,
    min_aspect_zone_area_cell_perc,
    include_buildings,
    outfolder,
):
    if h3_resolution:
        config["h3_resolution"] = h3_resolution

    if elev_diff_thresh:
        config["elev_diff_thresh"] = elev_diff_thresh

    if cell_buffer:
        config["cell_buffer"] = cell_buffer

    if min_area_urban:
        config["min_area_urban"] = min_area_urban

    if min_area_green:
        config["min_area_green"] = min_area_green

    if largest_urban:
        config["largest_urban"] = largest_urban

    if largest_green:
        config["largest_green"] = largest_green

    if min_aspect_zone_area_cell_perc:
        config["min_aspect_zone_area_cell_perc"] = min_aspect_zone_area_cell_perc

    if include_buildings:
        config["include_buildings"] = include_buildings

    if outfolder:
        config["outfolder"] = outfolder

    return config


def categoricals_to_str(
    df: gpd.GeoDataFrame,
) -> gpd.GeoDataFrame:
    for c in df.select_dtypes(include="category"):
        df[c] = df[c].astype(str)  # type:ignore

    df = gpd.GeoDataFrame(df)

    return df
