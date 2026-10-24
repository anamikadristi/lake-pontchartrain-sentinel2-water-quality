"""Step 3: score every candidate acquisition by the clear (cloud-free) fraction over
open water of each lake, using only the Scene Classification Layer (SCL) read at
160 m from the COG overviews. Output: data/s2_scene_scores.csv"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
from concurrent.futures import ThreadPoolExecutor
import numpy as np
import pandas as pd
from rasterio.enums import Resampling
from config import DATA, SCORE_BOXES, SCL_BAD
from lib_s2 import load_items, read_band, lonlat_to_rc

RES = 160.0
cat = pd.read_csv(DATA / "s2_catalog_15RYP.csv")
cand = cat[cat.cloud_cover <= 85].copy()
items = load_items()
boxes = {}
for k, (w, s, e, n) in SCORE_BOXES.items():
    r1, c0 = lonlat_to_rc(w, s, RES)
    r0, c1 = lonlat_to_rc(e, n, RES)
    boxes[k] = (slice(r0, r1), slice(c0, c1))


def score(iid):
    try:
        scl = read_band(items[iid], "SCL", RES, Resampling.nearest)
    except Exception as ex:  # network hiccup: record and move on
        return dict(id=iid, error=str(ex)[:120])
    out = dict(id=iid)
    for k, (rs, cs) in boxes.items():
        sub = scl[rs, cs]
        out[f"clear_{k}"] = float(np.isin(sub, list(SCL_BAD), invert=True).mean())
        out[f"sclwater_{k}"] = float((sub == 6).mean())
    return out


with ThreadPoolExecutor(8) as ex:
    res = list(ex.map(score, cand.id))
sc = cand.merge(pd.DataFrame(res), on="id", how="left")
sc.to_csv(DATA / "s2_scene_scores.csv", index=False)
print(len(cat), "acquisitions,", len(cand), "scored;",
      "errors:", sc["error"].notna().sum() if "error" in sc else 0)
print(sc[["clear_Pontchartrain", "clear_Maurepas"]].describe())
