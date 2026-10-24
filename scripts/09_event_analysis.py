"""Step 9: event analysis.
Plume definition: a pixel is "high turbidity" if its turbidity proxy exceeds the 90th
percentile of all Pontchartrain open-water pixels in the reference monthly scenes
(2021-2024, excluding 29 Aug - 31 Oct 2021 after Hurricane Ida; there were no spillway
openings 2021-2024). Plume fraction = high-turbidity pixels / valid pixels in the lake mask.
Also reported: lake-median turbidity, its ratio to the same-month reference median, and the
median in the "outfall sector" (open water within 15 km of the spillway outfall).
Outputs: outputs/tables/event_scene_metrics.csv, event_period_summary.csv,
         event_tests.csv, fig04_event_timeseries.png"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import json
import numpy as np
import pandas as pd
import rasterio
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.dates as mdates
from config import DATA, WORK, TAB, FIG, RES, BCS_OUTFALL, IDA_LANDFALL, BONNET_CARRE_2019
from lib_s2 import lonlat_to_rc, scene_indices

with rasterio.open(DATA / "lake_mask_40m.tif") as src:
    mask = src.read(1)
area = pd.read_csv(DATA / "lake_mask_areas.csv").set_index("lake").analysis_area_km2
st = pd.read_csv(TAB / "lake_scene_stats.csv", parse_dates=["date"])
u = st[st.used].copy()
r0, c0 = lonlat_to_rc(*BCS_OUTFALL, RES)
yy, xx = np.indices(mask.shape)
sector = (mask == 1) & (np.hypot(yy - r0, xx - c0) * RES <= 15000)
print("outfall sector open water km2:", round(sector.sum() * RES * RES / 1e6, 1))


def load(iid):
    return scene_indices(iid, ("turb",))["turb"]


ida_lo, ida_hi = pd.Timestamp(IDA_LANDFALL), pd.Timestamp("2021-10-31")
def is_ref(d):
    return (d.year >= 2021) and not (ida_lo <= d <= ida_hi)
P = u[u.lake == "Pontchartrain"]
refP = P[P.monthly & P.date.map(is_ref)]
pool = np.concatenate([(lambda a: a[(mask == 1) & np.isfinite(a)])(load(i)) for i in refP.id])
T90 = float(np.percentile(pool, 90))
print(f"reference scenes: {len(refP)}, pixels pooled: {pool.size}, P90 threshold = {T90:.2f}")

clim = pd.read_csv(TAB / "lake_monthly_climatology_2021_2024.csv")
cl = {(r.lake, r.month): r.turb_median for r in clim.itertuples()}
rows = []
for r in u.itertuples():
    a = load(r.id)
    code = 1 if r.lake == "Pontchartrain" else 2
    v = a[(mask == code) & np.isfinite(a)]
    d = dict(id=r.id, date=r.date, lake=r.lake, monthly=r.monthly, valid_frac=r.valid_frac,
             turb_median=r.turb_median, plume_frac=float((v > T90).mean()))
    d["plume_area_km2"] = d["plume_frac"] * area[r.lake]
    d["ratio_to_month_ref"] = r.turb_median / cl.get((r.lake, r.date.month), np.nan)
    if code == 1:
        s = a[sector & np.isfinite(a)]
        rest = a[(mask == 1) & ~sector & np.isfinite(a)]
        good = s.size > 0.5 * sector.sum() and rest.size > 0
        d["outfall_sector_median"] = float(np.median(s)) if good else np.nan
        # spatial contrast: outfall sector vs the rest of the lake on the same day,
        # which removes most of the lake-wide seasonal and wind signal
        d["sector_to_rest_ratio"] = float(np.median(s) / np.median(rest)) if good else np.nan
    rows.append(d)
em = pd.DataFrame(rows).sort_values(["lake", "date"])
em.to_csv(TAB / "event_scene_metrics.csv", index=False, float_format="%.4f")

periods = [("BCS 2019: before (1 Jan - 26 Feb 2019)", "2019-01-01", "2019-02-26"),
           ("BCS 2019: opening 1 (27 Feb - 11 Apr)", "2019-02-27", "2019-04-11"),
           ("BCS 2019: between openings (12 Apr - 9 May)", "2019-04-12", "2019-05-09"),
           ("BCS 2019: opening 2 (10 May - 27 Jul)", "2019-05-10", "2019-07-27"),
           ("BCS 2019: after (28 Jul - 30 Sep)", "2019-07-28", "2019-09-30"),
           ("Ida 2021: before (1 Jul - 28 Aug)", "2021-07-01", "2021-08-28"),
           ("Ida 2021: 1-14 days after (30 Aug - 12 Sep)", "2021-08-30", "2021-09-12"),
           ("Ida 2021: 15-63 days after (13 Sep - 31 Oct)", "2021-09-13", "2021-10-31")]
out = []
for name, a, b in periods:
    for lake, g in em[(em.date >= a) & (em.date <= b)].groupby("lake"):
        out.append(dict(period=name, lake=lake, n_scenes=len(g),
                        dates=" ".join(g.date.dt.strftime("%m-%d")),
                        turb_median=g.turb_median.median(),
                        ratio_to_month_ref=g.ratio_to_month_ref.median(),
                        plume_frac=g.plume_frac.median(), plume_area_km2=g.plume_area_km2.median(),
                        outfall_sector_median=g.outfall_sector_median.median() if lake == "Pontchartrain" else np.nan,
                        sector_to_rest_ratio=g.sector_to_rest_ratio.median() if lake == "Pontchartrain" else np.nan))
ps = pd.DataFrame(out)
ps.to_csv(TAB / "event_period_summary.csv", index=False, float_format="%.3f")
pd.set_option("display.width", 250)
print(ps.drop(columns="dates").to_string())

# Tests: Pontchartrain, 2019 opening months vs reference scenes of the same months
tests = []
Pm = em[em.lake == "Pontchartrain"]
open19 = Pm[Pm.date.apply(lambda d: any(pd.Timestamp(x) <= d <= pd.Timestamp(y) for x, y in BONNET_CARRE_2019))]
months = sorted(open19.date.dt.month.unique())
refm = Pm[Pm.date.map(is_ref) & Pm.monthly & Pm.date.dt.month.isin(months)]
for col in ["turb_median", "ratio_to_month_ref", "plume_frac", "outfall_sector_median", "sector_to_rest_ratio"]:
    x, y = open19[col].dropna(), refm[col].dropna()
    U, p = stats.mannwhitneyu(x, y, alternative="two-sided")
    tests.append(dict(test=f"Pontchartrain {col}: 2019 spillway-open scenes vs 2021-2024 same-month reference",
                      n_event=len(x), n_ref=len(y), median_event=x.median(), median_ref=y.median(),
                      mannwhitney_U=U, p_value=p))
# NDCI (chlorophyll-a proxy) after the 2019 openings: Jun-Aug 2019 vs Jun-Aug 2021-2024
for lake in ["Pontchartrain", "Maurepas"]:
    g = u[u.lake == lake]
    summer = g.date.dt.month.isin([6, 7, 8])
    ev = g[summer & (g.date.dt.year == 2019)].ndci_median
    rf = g[summer & g.date.map(is_ref) & g.monthly].ndci_median
    if len(ev) >= 2 and len(rf) >= 2:
        U, p = stats.mannwhitneyu(ev, rf, alternative="two-sided")
        tests.append(dict(test=f"{lake} lake-median NDCI: Jun-Aug 2019 scenes vs Jun-Aug 2021-2024 monthly scenes",
                          n_event=len(ev), n_ref=len(rf), median_event=ev.median(), median_ref=rf.median(),
                          mannwhitney_U=U, p_value=p))
for lake in ["Pontchartrain", "Maurepas"]:
    g = em[em.lake == lake]
    pre = g[(g.date >= "2021-07-01") & (g.date <= "2021-08-28")].turb_median
    post = g[(g.date >= "2021-08-30") & (g.date <= "2021-09-30")].turb_median
    if len(pre) >= 2 and len(post) >= 2:
        U, p = stats.mannwhitneyu(post, pre, alternative="two-sided")
    else:
        U, p = np.nan, np.nan
    tests.append(dict(test=f"{lake} lake-median turbidity: Ida post (30 Aug-30 Sep) vs pre (1 Jul-28 Aug 2021)",
                      n_event=len(post), n_ref=len(pre), median_event=post.median(),
                      median_ref=pre.median(), mannwhitney_U=U, p_value=p))
tt = pd.DataFrame(tests)
tt.to_csv(TAB / "event_tests.csv", index=False, float_format="%.4g")
print(tt.to_string())
json.dump(dict(T90=T90, n_ref_scenes=int(len(refP)), n_pixels=int(pool.size)),
          open(TAB / "plume_threshold.json", "w"), indent=1)

fig, axes = plt.subplots(1, 2, figsize=(14, 4.6))
for ax, (a, b), ttl in [(axes[0], ("2019-01-01", "2019-09-30"), "2019: Bonnet Carre Spillway openings"),
                        (axes[1], ("2021-07-01", "2021-10-31"), "2021: Hurricane Ida")]:
    for lake, c in [("Pontchartrain", "#1f5aa6"), ("Maurepas", "#c46a1b")]:
        g = em[(em.lake == lake) & (em.date >= a) & (em.date <= b)]
        ax.plot(g.date, g.plume_frac * 100, "o-", ms=4, lw=1, color=c, label=f"{lake}: % of open water > P90")
    if ttl.startswith("2019"):
        for x, y in BONNET_CARRE_2019:
            ax.axvspan(pd.Timestamp(x), pd.Timestamp(y), color="tan", alpha=0.4, lw=0)
        g = em[(em.lake == "Pontchartrain") & (em.date >= a) & (em.date <= b)]
        ax2 = ax.twinx()
        ax2.plot(g.date, g.sector_to_rest_ratio, "s--", ms=3.5, lw=0.8, color="0.35")
        ax2.axhline(np.nanmedian(em[(em.lake == "Pontchartrain") & em.date.map(is_ref) & em.monthly].sector_to_rest_ratio),
                    color="0.35", lw=0.6, ls=":")
        ax2.set_ylabel("outfall sector / rest-of-lake turbidity ratio\n(grey squares; dotted = median of all 2021-2024 reference scenes)", fontsize=8)
    else:
        ax.axvline(pd.Timestamp(IDA_LANDFALL), color="crimson", ls="--", lw=1)
    ax.set_title(ttl, fontsize=10)
    ax.set_ylabel(f"% of valid open water with turbidity proxy > {T90:.1f} (P90 of 2021-2024 reference)", fontsize=8)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8, loc="upper left")
    ax.xaxis.set_major_locator(mdates.MonthLocator())
    ax.xaxis.set_major_formatter(mdates.DateFormatter("%b %Y"))
    ax.tick_params(axis="x", labelsize=7, rotation=30)
plt.tight_layout()
plt.savefig(FIG / "fig04_event_plume_fraction.png", dpi=130)
