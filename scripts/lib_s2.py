"""Helpers: open Planetary Computer COGs and read AOI windows at reduced resolution."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import time
import numpy as np
import pystac
import planetary_computer
import rasterio
from rasterio.enums import Resampling
from rasterio.windows import from_bounds
from rasterio.transform import from_origin
from pyproj import Transformer
from config import AOI_LONLAT, EPSG, WORK, SCL_BAD

BANDS = ["B02", "B03", "B04", "B05", "B8A", "B11"]

_tr = Transformer.from_crs(4326, EPSG, always_xy=True)


def aoi_utm(res):
    w, s, e, n = AOI_LONLAT
    xs, ys = zip(*[_tr.transform(x, y) for x, y in [(w, s), (w, n), (e, s), (e, n)]])
    x0, x1 = np.floor(min(xs) / res) * res, np.ceil(max(xs) / res) * res
    y0, y1 = np.floor(min(ys) / res) * res, np.ceil(max(ys) / res) * res
    return x0, y0, x1, y1


def grid(res):
    x0, y0, x1, y1 = aoi_utm(res)
    width, height = int(round((x1 - x0) / res)), int(round((y1 - y0) / res))
    return from_origin(x0, y1, res, res), width, height


def load_items():
    ic = pystac.ItemCollection.from_file(str(WORK / "s2_items_15RYP.json"))
    return {it.id: it for it in ic}


def read_band(item, band, res, resampling=Resampling.average, tries=4):
    for k in range(tries):
        try:
            return _read_band(item, band, res, resampling)
        except Exception:
            if k == tries - 1:
                raise
            time.sleep(5 * (k + 1))


def _read_band(item, band, res, resampling=Resampling.average):
    """Stream only the AOI window of one band from its COG, resampled to `res` metres.
    GDAL uses the internal overviews, so only a small part of the file is fetched."""
    # drop any stale SAS token saved with the item, then sign with a fresh one
    href = planetary_computer.sign(item.assets[band].href.split("?")[0])
    transform, width, height = grid(res)
    x0, y0, x1, y1 = aoi_utm(res)
    with rasterio.Env(GDAL_DISABLE_READDIR_ON_OPEN="EMPTY_DIR", GDAL_HTTP_MAX_RETRY="4",
                      GDAL_HTTP_RETRY_DELAY="2", GDAL_HTTP_TIMEOUT="60",
                      GDAL_HTTP_CONNECTTIMEOUT="30", VSI_CACHE="TRUE"):
        with rasterio.open(href) as src:
            assert src.crs.to_epsg() == EPSG
            win = from_bounds(x0, y0, x1, y1, src.transform)
            arr = src.read(1, window=win, out_shape=(height, width), boundless=True,
                           fill_value=0, resampling=resampling)
    return arr


def lonlat_to_rc(lon, lat, res):
    transform, width, height = grid(res)
    x, y = _tr.transform(lon, lat)
    col = int((x - transform.c) / res)
    row = int((transform.f - y) / res)
    return row, col


def dogliotti(red, nir):
    """Dogliotti et al. (2015) turbidity (FNU) from water reflectance, switching red/NIR."""
    with np.errstate(divide="ignore", invalid="ignore"):
        t_red = 228.1 * red / (1 - red / 0.1641)
        t_nir = 3078.9 * nir / (1 - nir / 0.2112)
    w = np.clip((red - 0.05) / (0.07 - 0.05), 0, 1)
    t = (1 - w) * t_red + w * t_nir
    t[~((red > 0) & (red < 0.16) & (nir < 0.21))] = np.nan
    return t


SWIR_MAX = 0.03   # water pixels with B11 above this are treated as glint/haze contaminated


def scene_indices(iid, keys=("turb", "ndci", "cdom", "mndwi", "rgb", "scl", "b11")):
    """Water-quality proxies for one processed scene.
    1. valid = SCL not in {0,1,3,8,9,10,11}; water = valid and MNDWI > 0
    2. residual sun-glint / aerosol correction: subtract B11 (1610 nm) from every band,
       because water-leaving reflectance at 1610 nm is close to zero even in turbid water;
       pixels with B11 > SWIR_MAX are dropped
    3. turbidity = Dogliotti (B04 red, B8A NIR); NDCI = (B05-B04)/(B05+B04);
       CDOM proxy = B02/B03."""
    f = np.load(WORK / "refl" / f"{iid}.npz")
    r = {b: f[b].astype("float32") for b in BANDS}
    scl = f["scl"]
    out = {}
    with np.errstate(divide="ignore", invalid="ignore"):
        mndwi = (r["B03"] - r["B11"]) / (r["B03"] + r["B11"])
        nodata = f["rgb"].max(axis=2) == 0
        valid = ~nodata & ~np.isin(scl, list(SCL_BAD))
        mndwi = np.where(np.isfinite(mndwi), mndwi, -1.0)  # not stored = clearly land
        water = valid & (mndwi > 0) & (r["B11"] <= SWIR_MAX)
        c = {b: r[b] - r["B11"] for b in BANDS[:-1]}
        if "turb" in keys:
            out["turb"] = dogliotti(c["B04"], c["B8A"])
        if "ndci" in keys:
            out["ndci"] = (c["B05"] - c["B04"]) / (c["B05"] + c["B04"])
        if "cdom" in keys:
            out["cdom"] = c["B02"] / c["B03"]
            out["cdom"][(c["B02"] <= 0) | (c["B03"] <= 0)] = np.nan
    for k in ("turb", "ndci", "cdom"):
        if k in out:
            out[k][~water] = np.nan
    if "ndci" in out:
        out["ndci"][(c["B04"] <= 0) | (c["B05"] <= 0)] = np.nan
    if "mndwi" in keys:
        mndwi[~valid] = np.nan
        out["mndwi"] = mndwi
    if "b11" in keys:
        b = r["B11"].copy()
        b[~(valid & (mndwi > 0))] = np.nan
        out["b11"] = b
    if "rgb" in keys:
        out["rgb"] = f["rgb"]
    if "scl" in keys:
        out["scl"] = scl
    return out
