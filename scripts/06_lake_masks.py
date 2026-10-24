"""Step 6: build a static open-water mask for each lake from the 72 monthly scenes.
  water frequency = (# clear obs with MNDWI > 0) / (# clear obs)
  open water      = frequency >= 0.9 and >= 12 clear observations
  shoreline buffer: erode 3 pixels (120 m) to reduce mixed land/water pixels and
                    adjacency effects; this also cuts the narrow passes, so each lake
                    becomes its own connected component, picked with a seed point.
Outputs: data/lake_mask_40m.tif (1 = Pontchartrain, 2 = Maurepas),
         data/lake_open_water.geojson (polygons, EPSG:4326), water-frequency figure."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import json
import numpy as np
import pandas as pd
import rasterio
from rasterio import features
from rasterio.warp import transform_geom
from scipy import ndimage
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from config import DATA, WORK, RES, EPSG, LAKE_SEEDS, FIG
from lib_s2 import grid, lonlat_to_rc, scene_indices

sel = pd.read_csv(DATA / "s2_selected_scenes.csv")
mon = sel[sel.monthly]
transform, W, H = grid(RES)
nwat = np.zeros((H, W), "int16")
nobs = np.zeros((H, W), "int16")
for iid in mon.id:
    m = scene_indices(iid, ("mndwi",))["mndwi"]
    ok = np.isfinite(m)
    nobs += ok
    nwat += ok & (m > 0)
freq = np.where(nobs > 0, nwat / np.maximum(nobs, 1), np.nan)
open_w = (freq >= 0.9) & (nobs >= 12)
open_w = ndimage.binary_opening(open_w, iterations=1)
core = ndimage.binary_erosion(open_w, iterations=3)
lab, n = ndimage.label(core)
mask = np.zeros((H, W), "uint8")
areas = {}
for code, (name, seeds) in enumerate(LAKE_SEEDS.items(), start=1):
    for lon, lat in seeds:
        r, c = lonlat_to_rc(lon, lat, RES)
        k = lab[r, c]
        assert k > 0, f"seed for {name} not in water"
        assert mask[lab == k].max() in (0, code), "lakes are still connected"
        mask[lab == k] = code
        print(name, (lon, lat), "component km2:", round((lab == k).sum() * RES * RES / 1e6, 1))
    areas[name] = float((mask == code).sum() * RES * RES / 1e6)
print("analysis area (km2):", {k: round(v, 1) for k, v in areas.items()})

prof = dict(driver="GTiff", height=H, width=W, count=1, dtype="uint8", crs=f"EPSG:{EPSG}",
            transform=transform, compress="deflate", nodata=0)
with rasterio.open(DATA / "lake_mask_40m.tif", "w", **prof) as dst:
    dst.write(mask, 1)
feats = []
for geom, val in features.shapes(mask, mask=mask > 0, transform=transform):
    g = transform_geom(f"EPSG:{EPSG}", "EPSG:4326", geom)
    name = list(LAKE_SEEDS)[int(val) - 1]
    feats.append(dict(type="Feature", geometry=g, properties=dict(lake=name, code=int(val))))
(DATA / "lake_open_water.geojson").write_text(json.dumps(dict(type="FeatureCollection", features=feats)))
pd.DataFrame([dict(lake=k, analysis_area_km2=round(v, 1)) for k, v in areas.items()]).to_csv(
    DATA / "lake_mask_areas.csv", index=False)

fig, ax = plt.subplots(figsize=(10, 5.2))
im = ax.imshow(freq, cmap="Blues", vmin=0, vmax=1)
ax.contour(mask > 0, levels=[0.5], colors="crimson", linewidths=0.6)
ax.set_title(f"Water frequency (MNDWI > 0) in {len(mon)} monthly Sentinel-2 scenes, 2019-2024;\n"
             "red outline: open-water analysis mask (120 m shoreline buffer)")
ax.set_axis_off()
plt.colorbar(im, ax=ax, shrink=0.7, label="fraction of clear observations")
plt.tight_layout()
plt.savefig(FIG / "fig01_water_frequency_mask.png", dpi=130)
