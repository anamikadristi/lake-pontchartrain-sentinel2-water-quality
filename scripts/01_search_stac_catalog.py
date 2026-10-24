"""Step 1: search the anonymous Microsoft Planetary Computer STAC for Sentinel-2 L2A
items over tile 15RYP (Lake Pontchartrain and Lake Maurepas), 2019-2024.
Writes data/s2_catalog_15RYP.csv (metadata only, no imagery)."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
import pystac_client
import planetary_computer
from config import STAC_URL, COLLECTION, MGRS_TILE, START, END, AOI_LONLAT, DATA, WORK

cat = pystac_client.Client.open(STAC_URL, modifier=planetary_computer.sign_inplace)
search = cat.search(collections=[COLLECTION], bbox=list(AOI_LONLAT),
                    datetime=f"{START}/{END}",
                    query={"s2:mgrs_tile": {"eq": MGRS_TILE}})
items = search.item_collection()
items.save_object(str(WORK / "s2_items_15RYP.json"))  # full STAC metadata, re-signed on use
rows = []
for it in items:
    p = it.properties
    rows.append(dict(id=it.id, datetime=p["datetime"], platform=p.get("platform"),
                     rel_orbit=p.get("sat:relative_orbit"),
                     cloud_cover=p.get("eo:cloud_cover"),
                     baseline=p.get("s2:processing_baseline"),
                     west=it.bbox[0], south=it.bbox[1], east=it.bbox[2], north=it.bbox[3]))
df = pd.DataFrame(rows).sort_values("datetime")
df["date"] = pd.to_datetime(df["datetime"]).dt.strftime("%Y-%m-%d")
n0 = len(df)
# keep one item per acquisition (the most recent processing if duplicated)
df = df.sort_values(["date", "rel_orbit", "id"]).drop_duplicates(["date", "rel_orbit"], keep="last")
df.to_csv(DATA / "s2_catalog_15RYP.csv", index=False)
print(f"{n0} items found, {len(df)} unique acquisitions, {df.date.min()} to {df.date.max()}")
print(df.groupby("rel_orbit").size())
