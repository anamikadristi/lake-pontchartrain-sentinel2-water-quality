"""Shared settings for the Lake Pontchartrain / Lake Maurepas Sentinel-2 workflow."""
import os
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
DATA = REPO / "data"
FIG = REPO / "outputs" / "figures"
TAB = REPO / "outputs" / "tables"
# Large intermediate rasters (not committed). Override with env var S2WQ_WORK.
WORK = Path(os.environ.get("S2WQ_WORK", REPO / "work"))
for p in (DATA, FIG, TAB, WORK):
    p.mkdir(parents=True, exist_ok=True)

STAC_URL = "https://planetarycomputer.microsoft.com/api/stac/v1"
COLLECTION = "sentinel-2-l2a"
MGRS_TILE = "15RYP"          # one tile covers both lakes (east edge about -89.77)
EPSG = 32615                  # UTM 15N, native CRS of tile 15RYP
START, END = "2019-01-01", "2024-12-31"

# Analysis window (lon/lat) and output grid
AOI_LONLAT = (-90.65, 29.99, -89.77, 30.40)
RES = 40.0                    # metres; read from COG overviews (10 m bands -> 40 m)

# Rough rectangles used only to score cloud cover over open water
SCORE_BOXES = {
    "Pontchartrain": (-90.38, 30.05, -89.82, 30.30),
    "Maurepas": (-90.57, 30.21, -90.43, 30.31),
}
# Seed points used to pick each lake from the connected water components. The Causeway
# and the I-10 twin span split Lake Pontchartrain into three pieces at 40 m.
LAKE_SEEDS = {"Pontchartrain": [(-90.25, 30.18), (-90.00, 30.18), (-89.80, 30.17)],
              "Maurepas": [(-90.50, 30.26)]}

# Event dates
BONNET_CARRE_2019 = [("2019-02-27", "2019-04-11"), ("2019-05-10", "2019-07-27")]
IDA_LANDFALL = "2021-08-29"
# Bonnet Carre Spillway outfall into Lake Pontchartrain (approx. lake-side end)
BCS_OUTFALL = (-90.43, 30.06)

SCL_BAD = {0, 1, 3, 8, 9, 10, 11}   # no data, saturated, shadow, cloud med/high, cirrus, snow
