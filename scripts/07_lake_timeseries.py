"""Step 7: lake-wide statistics for every processed scene inside the static open-water
masks. A lake/scene pair is kept only if >= 50% of the mask has valid (cloud-free,
MNDWI > 0) pixels, so partly cloudy scenes do not bias the lake median.
Outputs: outputs/tables/lake_scene_stats.csv, lake_monthly_climatology.csv, fig02."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
import rasterio
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from config import DATA, WORK, TAB, FIG, BONNET_CARRE_2019, IDA_LANDFALL
from lib_s2 import scene_indices

sel = pd.read_csv(DATA / "s2_selected_scenes.csv", parse_dates=["date"])
with rasterio.open(DATA / "lake_mask_40m.tif") as src:
    mask = src.read(1)
LAKES = {1: "Pontchartrain", 2: "Maurepas"}
rows = []
for r in sel.itertuples():
    arr = scene_indices(r.id, ("turb", "ndci", "cdom", "b11"))
    b11 = arr.pop("b11")
    for code, lake in LAKES.items():
        m = mask == code
        t = arr["turb"][m]
        ok = np.isfinite(t)
        d = dict(id=r.id, date=r.date, monthly=r.monthly, reason=r.reason, lake=lake,
                 valid_frac=ok.mean(), b11_median=np.nanmedian(b11[m]) if np.isfinite(b11[m]).any() else np.nan)
        for k, a in arr.items():
            v = a[m][ok]
            v = v[np.isfinite(v)]
            if v.size:
                d[f"{k}_median"] = np.median(v)
                d[f"{k}_p10"], d[f"{k}_p90"] = np.percentile(v, [10, 90])
                d[f"{k}_mean"] = v.mean()
        rows.append(d)
st = pd.DataFrame(rows)
st["used"] = st.valid_frac >= 0.5
st.to_csv(TAB / "lake_scene_stats.csv", index=False, float_format="%.4f")
u = st[st.used]
print("lake-scene pairs used:", u.groupby("lake").size().to_dict())

# Reference climatology: 2021-2024, no spillway openings (openings: 2019 and Apr 3-May 1 2020),
# excluding Hurricane Ida (Aug 29 - Oct 31 2021)
ref = u[(u.date.dt.year >= 2021) & ~((u.date >= IDA_LANDFALL) & (u.date <= "2021-10-31"))].copy()
ref["month"] = ref.date.dt.month
clim = ref.groupby(["lake", "month"]).agg(n=("turb_median", "size"),
                                          turb_median=("turb_median", "median"),
                                          ndci_median=("ndci_median", "median"),
                                          cdom_median=("cdom_median", "median")).reset_index()
clim.to_csv(TAB / "lake_monthly_climatology_2021_2024.csv", index=False, float_format="%.4f")
summ = u.groupby("lake").agg(n_scenes=("id", "size"), first=("date", "min"), last=("date", "max"),
                             turb_median_of_medians=("turb_median", "median"),
                             turb_min=("turb_median", "min"), turb_max=("turb_median", "max"),
                             ndci_median=("ndci_median", "median"), cdom_median=("cdom_median", "median"))
summ.to_csv(TAB / "lake_summary.csv", float_format="%.4f")
print(summ.to_string())

fig, axes = plt.subplots(3, 1, figsize=(12, 9), sharex=True)
labels = {"turb": "Turbidity proxy (FNU-like,\nDogliotti 2015, uncalibrated)",
          "ndci": "NDCI (chlorophyll-a proxy)", "cdom": "B2/B3 ratio (CDOM proxy;\nlower = more absorption)"}
colors = {"Pontchartrain": "#1f5aa6", "Maurepas": "#c46a1b"}
for ax, k in zip(axes, labels):
    for lake, g in u.groupby("lake"):
        g = g.sort_values("date")
        ax.plot(g.date, g[f"{k}_median"], "-", color=colors[lake], lw=0.8, alpha=0.6)
        ax.plot(g.date, g[f"{k}_median"], "o", ms=2.5, color=colors[lake], label=lake)
    for a, b in BONNET_CARRE_2019 + [("2020-04-03", "2020-05-01")]:
        ax.axvspan(pd.Timestamp(a), pd.Timestamp(b), color="tan", alpha=0.35, lw=0)
    ax.axvline(pd.Timestamp(IDA_LANDFALL), color="crimson", lw=1, ls="--")
    ax.set_ylabel(labels[k], fontsize=9)
    ax.grid(alpha=0.3)
axes[0].set_yscale("log")
axes[0].legend(loc="upper right", fontsize=8)
axes[0].set_title("Lake-median Sentinel-2 water-quality proxies, 2019-2024 (shaded: Bonnet Carre "
                  "Spillway open; dashed: Hurricane Ida landfall)", fontsize=10)
axes[-1].xaxis.set_major_locator(mdates.YearLocator())
plt.tight_layout()
plt.savefig(FIG / "fig02_lake_timeseries.png", dpi=130)
