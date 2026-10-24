"""Step 5: for each selected acquisition, stream only the AOI window of bands B02, B03,
B04, B05, B8A, B11 and SCL from the Planetary Computer COGs, resampled to 40 m (GDAL
reads the internal overviews; no full tile is downloaded).
Reflectance = (DN + BOA_ADD_OFFSET) / 10000, with offset -1000 for processing baseline
>= 04.00 (25 Jan 2022 onward). Reflectances are stored as float16 for pixels that are
not clearly land (MNDWI > -0.2), plus SCL and an 8-bit true-colour quicklook.
Indices are computed later by lib_s2.scene_indices(), so the algorithms can be changed
without downloading again. Output: WORK/refl/<id>.npz"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import os
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
import numpy as np
import pandas as pd
from rasterio.enums import Resampling
from config import DATA, WORK, RES
from lib_s2 import load_items, read_band, BANDS

OUT = WORK / "refl"
OUT.mkdir(parents=True, exist_ok=True)


def process(row):
    out = OUT / f"{row.id}.npz"
    if out.exists():
        return row.id, "cached"
    item = ITEMS[row.id]
    t0 = time.time()
    with ThreadPoolExecutor(7) as ex:
        futs = {b: ex.submit(read_band, item, b, RES) for b in BANDS}
        fscl = ex.submit(read_band, item, "SCL", RES, Resampling.nearest)
        dn = {b: f.result().astype("float32") for b, f in futs.items()}
        scl = fscl.result().astype("uint8")
    off = -1000.0 if float(row.baseline) >= 4.0 else 0.0
    nodata = dn["B03"] == 0
    rho = {b: (v + off) / 10000.0 for b, v in dn.items()}
    with np.errstate(divide="ignore", invalid="ignore"):
        mndwi = (rho["B03"] - rho["B11"]) / (rho["B03"] + rho["B11"])
    keep = ~nodata & (mndwi > -0.2)
    rgb = np.dstack([rho["B04"], rho["B03"], rho["B02"]])
    rgb = (np.clip(rgb / 0.18, 0, 1) ** (1 / 1.8) * 255).astype("uint8")
    rgb[nodata] = 0
    arrs = {b: np.where(keep, np.clip(v, -0.5, 1.5), np.nan).astype("float16") for b, v in rho.items()}
    tmp = OUT / f"{row.id}.{os.getpid()}.tmp.npz"
    np.savez_compressed(tmp, scl=scl, rgb=rgb, offset=off, **arrs)
    tmp.replace(out)  # atomic: an interrupted run never leaves a partial file
    return row.id, f"{time.time() - t0:.0f}s"


if __name__ == "__main__":
    ITEMS = load_items()
    sel = pd.read_csv(DATA / "s2_selected_scenes.csv")
    if len(sys.argv) > 1 and sys.argv[1] == "--reverse":  # lets a 2nd process work from the end
        sel = sel.iloc[::-1]
    elif len(sys.argv) > 1:
        sel = sel.head(int(sys.argv[1]))
    with ThreadPoolExecutor(3) as ex:
        futs = [ex.submit(process, r) for r in sel.itertuples()]
        for i, f in enumerate(as_completed(futs)):
            try:
                iid, msg = f.result()
            except Exception as e:
                iid, msg = "?", f"FAILED {e}"
            print(i + 1, iid, msg, flush=True)
