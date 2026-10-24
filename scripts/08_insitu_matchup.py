"""Step 8: satellite vs in-situ matchups.
Satellite value: median of valid water pixels (cloud-free, MNDWI > 0) in the 3 x 3 window
(120 m x 120 m at 40 m) centred on the station; at least 5 valid pixels required.
  (a) LDEQ discrete turbidity (NTU, WQP), shallowest sample per station and date,
      scene within +/- 1 day (same-day subset also reported)
  (b) USGS 300138090064700 New Canal Lighthouse continuous turbidity (FNU), median of
      30-min values within +/- 1 h of the overpass
  (c) chlorophyll: NDCI vs USGS discrete chlorophyll-a (ug/L, +/- 1 day) and vs USGS
      continuous chlorophyll fluorescence (+/- 1 h)
Outputs: outputs/tables/matchups_*.csv, matchup_statistics.csv, fig03."""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import numpy as np
import pandas as pd
from scipy import stats
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from config import DATA, WORK, TAB, FIG, RES
from lib_s2 import load_items, lonlat_to_rc, scene_indices

sel = pd.read_csv(DATA / "s2_selected_scenes.csv", parse_dates=["date"])
items = load_items()
sel["time"] = [pd.Timestamp(items[i].properties["datetime"]) for i in sel.id]
_cache = {}


def scene(iid):
    if iid not in _cache:
        if len(_cache) > 6:
            _cache.clear()
        _cache[iid] = scene_indices(iid, ("turb", "ndci"))
    return _cache[iid]


def sat_value(iid, lon, lat, var):
    r, c = lonlat_to_rc(lon, lat, RES)
    a = scene(iid)[var]
    if not (1 <= r < a.shape[0] - 1 and 1 <= c < a.shape[1] - 1):
        return np.nan, 0
    w = a[r - 1:r + 2, c - 1:c + 2]
    w = w[np.isfinite(w)]
    return (float(np.median(w)) if w.size >= 5 else np.nan), int(w.size)


def match_discrete(df, var, max_days=1):
    out = []
    for o in df.itertuples():
        for s in sel[(sel.date - o.date).abs() <= pd.Timedelta(days=max_days)].itertuples():
            v, n = sat_value(s.id, o.lon, o.lat, var)
            if np.isfinite(v):
                out.append(dict(station=o.station, name=o.name, lon=o.lon, lat=o.lat,
                                insitu_date=o.date.date(), insitu=o.value, scene=s.id,
                                scene_date=s.date.date(), dt_days=(s.date - o.date).days,
                                sat=v, n_pix=n))
    m = pd.DataFrame(out)
    if len(m):  # one matchup per observation: closest scene in time
        m["adt"] = m.dt_days.abs()
        m = m.sort_values("adt").groupby(["station", "insitu_date"]).head(1).drop(columns="adt")
    return m


def match_continuous(cont, var, lon, lat):
    out = []
    for s in sel.itertuples():
        w = cont[(cont.time - s.time).abs() <= pd.Timedelta(hours=1)]
        if len(w) < 2:
            continue
        v, n = sat_value(s.id, lon, lat, var)
        if np.isfinite(v):
            out.append(dict(scene=s.id, scene_time=s.time, insitu=w.value.median(),
                            n_insitu=len(w), sat=v, n_pix=n))
    return pd.DataFrame(out)


def mstats(m, label, log=True):
    x, y = m.insitu.values, m.sat.values
    d = dict(comparison=label, n=len(m))
    if len(m) < 4:
        return d
    d["spearman_rho"], d["spearman_p"] = stats.spearmanr(x, y)
    if log:
        ok = (x > 0) & (y > 0)
        lx, ly = np.log10(x[ok]), np.log10(y[ok])
        d["pearson_r_log10"], d["pearson_p_log10"] = stats.pearsonr(lx, ly)
        lr = stats.linregress(lx, ly)
        d["slope_log10"], d["intercept_log10"] = lr.slope, lr.intercept
        d["rmse_log10"] = float(np.sqrt(np.mean((ly - lx) ** 2)))
        d["median_ratio_sat_over_insitu"] = float(np.median(y[ok] / x[ok]))
        d["median_abs_pct_diff"] = float(np.median(np.abs(y[ok] - x[ok]) / x[ok]) * 100)
        d["insitu_median"], d["sat_median"] = float(np.median(x)), float(np.median(y))
    else:
        d["pearson_r"], d["pearson_p"] = stats.pearsonr(x, y)
    return d


# ---- (a) LDEQ discrete turbidity
wq = pd.read_csv(DATA / "insitu_wqp_lake_stations.csv")
wq = wq.rename(columns={"MonitoringLocationIdentifier": "station", "MonitoringLocationName": "name",
                        "LatitudeMeasure": "lat", "LongitudeMeasure": "lon"})
wq["date"] = pd.to_datetime(wq.ActivityStartDate)
wq["depth"] = wq["ActivityDepthHeightMeasure/MeasureValue"].fillna(0)
tb = wq[wq.CharacteristicName == "Turbidity"].sort_values("depth")
tb = tb.groupby(["station", "date"]).head(1)
m_ldeq = match_discrete(tb, "turb")
m_ldeq.to_csv(TAB / "matchups_ldeq_turbidity.csv", index=False, float_format="%.3f")

# ---- (b) USGS continuous turbidity at New Canal Lighthouse
NC = (-90.1131667, 30.0271667)
ct = pd.read_csv(WORK / "insitu" / "usgs_continuous_63680.csv")
ct["time"] = pd.to_datetime(ct.time, utc=True)
ct = ct.dropna(subset=["value"])
m_nc = match_continuous(ct, "turb", *NC)
m_nc.to_csv(TAB / "matchups_usgs_newcanal_turbidity.csv", index=False, float_format="%.3f")

# ---- (c) chlorophyll
ch = wq[wq.CharacteristicName == "Chlorophyll a"].sort_values("depth").groupby(["station", "date"]).head(1)
m_chl = match_discrete(ch, "ndci")
m_chl.to_csv(TAB / "matchups_usgs_chla_ndci.csv", index=False, float_format="%.4f")
cf = pd.read_csv(WORK / "insitu" / "usgs_continuous_32316.csv")
cf["time"] = pd.to_datetime(cf.time, utc=True)
cf = cf.dropna(subset=["value"])
m_cf = match_continuous(cf, "ndci", *NC)
m_cf.to_csv(TAB / "matchups_usgs_newcanal_chlfluor_ndci.csv", index=False, float_format="%.4f")

rows = [mstats(m_ldeq, "Turbidity: S2 proxy vs LDEQ discrete (NTU), +/-1 day"),
        mstats(m_ldeq[m_ldeq.dt_days == 0], "Turbidity: S2 proxy vs LDEQ discrete (NTU), same day"),
        mstats(m_nc, "Turbidity: S2 proxy vs USGS New Canal continuous (FNU), +/-1 h"),
        mstats(m_chl, "NDCI vs USGS discrete chlorophyll-a (ug/L), +/-1 day", log=False),
        mstats(m_cf, "NDCI vs USGS New Canal chlorophyll fluorescence (ug/L), +/-1 h", log=False)]
if len(m_ldeq) and len(m_nc):
    allt = pd.concat([m_ldeq[["insitu", "sat"]], m_nc[["insitu", "sat"]]])
    rows.append(mstats(allt, "Turbidity: all in-situ pooled"))
ms = pd.DataFrame(rows)
ms.to_csv(TAB / "matchup_statistics.csv", index=False, float_format="%.3f")
pd.set_option("display.width", 250)
print(ms.T.to_string())

fig, axes = plt.subplots(1, 3, figsize=(14, 4.6))
for ax, m, ttl, c in [(axes[0], m_ldeq, "LDEQ discrete turbidity (NTU)", "#1f5aa6"),
                      (axes[1], m_nc, "USGS New Canal continuous turbidity (FNU)", "#2a8a4a")]:
    if len(m):
        ax.loglog(m.insitu, m.sat, "o", ms=5, color=c, alpha=0.75)
        lo = max(0.5, min(m.insitu.min(), m.sat.min()) * 0.7)
        hi = max(m.insitu.max(), m.sat.max()) * 1.4
        ax.plot([lo, hi], [lo, hi], "k--", lw=0.8, label="1:1")
        r = ms[ms.comparison.str.contains(ttl.split()[0]) & ms.comparison.str.contains("Turbidity")].iloc[0]
        ax.text(0.04, 0.96, f"n = {int(r.n)}\nSpearman rho = {r.spearman_rho:.2f}\n"
                f"median ratio S2/in situ = {r.median_ratio_sat_over_insitu:.2f}",
                transform=ax.transAxes, va="top", fontsize=8.5)
        ax.legend(loc="lower right", fontsize=8)
    ax.set_xlabel(ttl)
    ax.set_ylabel("Sentinel-2 turbidity proxy (FNU-like)")
    ax.grid(alpha=0.3, which="both")
ax = axes[2]
for m, lbl, mk in [(m_chl, "USGS discrete chl-a", "o"), (m_cf, "USGS New Canal chl fluorescence", "^")]:
    if len(m):
        ax.plot(m.insitu, m.sat, mk, ms=5, alpha=0.75, label=f"{lbl} (n={len(m)})")
ax.set_xlabel("in-situ chlorophyll (ug/L)")
ax.set_ylabel("NDCI (unitless proxy)")
ax.legend(fontsize=8)
ax.grid(alpha=0.3)
fig.suptitle("Sentinel-2 proxies vs public in-situ data (3 x 3 pixel median at 40 m)", fontsize=11)
plt.tight_layout()
plt.savefig(FIG / "fig03_insitu_matchups.png", dpi=130)
