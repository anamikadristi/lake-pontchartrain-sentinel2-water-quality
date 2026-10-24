"""Step 10: maps, animated GIF and GIS-ready rasters.
  fig05: pixel-wise median of the reference monthly scenes (2021-2024, Ida excluded)
         for the three proxies
  fig06: Bonnet Carre 2019 sequence (true colour and turbidity proxy)
  fig07: Hurricane Ida 2021 sequence
  outputs/figures/turbidity_monthly_2019_2024.gif: one frame per monthly scene, 120 m
  outputs/rasters/*.tif: GeoTIFF (EPSG:32615) for QGIS / ArcGIS Pro"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import io
import numpy as np
import pandas as pd
import rasterio
import imageio.v2 as imageio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LogNorm
from config import DATA, WORK, TAB, FIG, REPO, RES, EPSG, IDA_LANDFALL
from lib_s2 import grid, scene_indices

RAS = REPO / "outputs" / "rasters"
RAS.mkdir(parents=True, exist_ok=True)
with rasterio.open(DATA / "lake_mask_40m.tif") as src:
    mask = src.read(1)
transform, W, H = grid(RES)
em = pd.read_csv(TAB / "event_scene_metrics.csv", parse_dates=["date"])
st = pd.read_csv(TAB / "lake_scene_stats.csv", parse_dates=["date"])
TN = LogNorm(2, 150)
TCMAP = "YlOrBr"


def load(iid, keys=("turb", "rgb")):
    return scene_indices(iid, keys)


def write_tif(path, arr):
    prof = dict(driver="GTiff", height=H, width=W, count=1, dtype="float32", crs=f"EPSG:{EPSG}",
                transform=transform, compress="deflate", predictor=3, nodata=np.nan, tiled=True)
    with rasterio.open(path, "w", **prof) as dst:
        dst.write(arr.astype("float32"), 1)


# ---- reference medians
P = st[(st.lake == "Pontchartrain") & st.used & st.monthly]
ida = (P.date >= IDA_LANDFALL) & (P.date <= "2021-10-31")
ref = P[(P.date.dt.year >= 2021) & ~ida]
inm = mask > 0
med = {}
for k in ("turb", "ndci", "cdom"):
    stack = np.stack([load(i, (k,))[k][inm] for i in ref.id])
    m = np.full(mask.shape, np.nan, "float32")
    v = np.nanmedian(stack, axis=0)
    v[np.isfinite(stack).sum(0) < 10] = np.nan   # need >= 10 clear observations
    m[inm] = v
    med[k] = m
write_tif(RAS / "turbidity_proxy_median_2021_2024_40m.tif", med["turb"])
write_tif(RAS / "ndci_median_2021_2024_40m.tif", med["ndci"])
fig, axes = plt.subplots(3, 1, figsize=(10, 13))
specs = [("turb", "Turbidity proxy (FNU-like)", TCMAP, TN),
         ("ndci", "NDCI (chlorophyll-a proxy)", "viridis", None),
         ("cdom", "B2/B3 ratio (CDOM proxy, lower = more absorption)", "cividis_r", None)]
for ax, (k, lbl, cm, nrm) in zip(axes, specs):
    v = med[k]
    kw = dict(norm=nrm) if nrm else dict(vmin=np.nanpercentile(v, 2), vmax=np.nanpercentile(v, 98))
    ax.imshow(np.where(inm, np.nan, 1), cmap="Greys", vmin=0, vmax=3)
    im = ax.imshow(v, cmap=cm, **kw)
    plt.colorbar(im, ax=ax, shrink=0.75, label=lbl)
    ax.set_axis_off()
axes[0].set_title(f"Pixel-wise median of {len(ref)} reference monthly scenes (2021-2024, Ida excluded)\n"
                  "open-water analysis mask only; proxies, not calibrated concentrations", fontsize=10)
plt.tight_layout()
plt.savefig(FIG / "fig05_reference_median_maps.png", dpi=110)
plt.close()


def best(a, b):
    g = em[(em.lake == "Pontchartrain") & (em.date >= a) & (em.date <= b)]
    return g.sort_values("valid_frac").iloc[-1] if len(g) else None


def sequence(panels, fname, title):
    sel = [(lbl, best(a, b)) for lbl, a, b in panels]
    sel = [(l, r) for l, r in sel if r is not None]
    fig, axes = plt.subplots(2, len(sel), figsize=(4.0 * len(sel) + 1, 5.0), layout="constrained")
    for j, (lbl, r) in enumerate(sel):
        d = load(r.id)
        axes[0, j].imshow(d["rgb"])
        axes[0, j].set_title(f"{lbl}\n{r.date:%Y-%m-%d}", fontsize=9)
        axes[1, j].imshow(d["rgb"].mean(2), cmap="Greys_r", vmin=0, vmax=400)
        im = axes[1, j].imshow(d["turb"], cmap=TCMAP, norm=TN)
        axes[1, j].set_xlabel(f"Pontchartrain > P90: {r.plume_frac * 100:.0f}% of open water", fontsize=8)
        for ax in axes[:, j]:
            ax.set_xticks([]); ax.set_yticks([])
    fig.colorbar(im, ax=axes[:, -1].tolist(), shrink=0.6, label="turbidity proxy (FNU-like)")
    fig.suptitle(title, fontsize=11)
    plt.savefig(FIG / fname, dpi=110, bbox_inches="tight")
    plt.close()
    return sel


s6 = sequence([("Before opening", "2019-01-01", "2019-02-26"),
               ("Opening 1", "2019-02-27", "2019-04-11"),
               ("Between openings", "2019-04-12", "2019-05-09"),
               ("Opening 2", "2019-05-10", "2019-07-27"),
               ("After closure", "2019-07-28", "2019-09-30")],
              "fig06_bonnet_carre_2019_maps.png",
              "Bonnet Carre Spillway 2019: clearest Sentinel-2 scene in each period (top: true colour; bottom: turbidity proxy on water)")
s7 = sequence([("Before Ida", "2021-07-01", "2021-08-28"),
               ("1-14 days after", "2021-08-30", "2021-09-12"),
               ("15-35 days after", "2021-09-13", "2021-10-03"),
               ("36-63 days after", "2021-10-04", "2021-10-31")],
              "fig07_hurricane_ida_2021_maps.png",
              "Hurricane Ida (landfall 29 Aug 2021): clearest scene in each period")
r = dict(s6)["Opening 2"]
write_tif(RAS / f"turbidity_proxy_{r.date:%Y%m%d}_bonnet_carre_opening2_40m.tif", load(r.id, ("turb",))["turb"])
pd.DataFrame([dict(figure="fig06", panel=l, scene=x.id, date=x.date.date()) for l, x in s6] +
             [dict(figure="fig07", panel=l, scene=x.id, date=x.date.date()) for l, x in s7]).to_csv(
    TAB / "map_scene_list.csv", index=False)

# ---- GIF of monthly scenes (decimated to 120 m)
sel = pd.read_csv(DATA / "s2_selected_scenes.csv", parse_dates=["date"])
frames = []
for rr in sel[sel.monthly].itertuples():
    d = load(rr.id)
    t = d["turb"][::3, ::3]
    g = d["rgb"][::3, ::3].mean(2)
    fig, ax = plt.subplots(figsize=(6.4, 3.75), dpi=80)
    ax.imshow(g, cmap="Greys_r", vmin=0, vmax=400)
    im = ax.imshow(t, cmap=TCMAP, norm=TN)
    ax.set_axis_off()
    ax.set_title(f"Sentinel-2 turbidity proxy  {rr.date:%Y-%m-%d}", fontsize=10)
    plt.colorbar(im, ax=ax, shrink=0.7, label="FNU-like")
    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format="png")
    plt.close()
    frames.append(imageio.imread(buf.getvalue()))
imageio.mimsave(FIG / "turbidity_monthly_2019_2024.gif", frames, duration=0.6, loop=0)
print("GIF frames:", len(frames), "size MB:", round((FIG / "turbidity_monthly_2019_2024.gif").stat().st_size / 1e6, 2))
