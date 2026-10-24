"""Step 4: choose a manageable scene list.
  monthly : the clearest acquisition per calendar month (weighted 85% Pontchartrain,
            15% Maurepas, roughly their open-water areas), if Pontchartrain is >= 50% clear
  event   : every acquisition with >= 50% clear Pontchartrain in Jan-Sep 2019
            (Bonnet Carre openings) and Jul-Oct 2021 (Hurricane Ida)
  matchup : acquisitions within +/- 1 day of a discrete in-situ turbidity or chlorophyll
            sample, or during the USGS New Canal continuous record, with >= 30% clear
            Pontchartrain (station-level cloud screening is done later per pixel)
Output: data/s2_selected_scenes.csv"""
import sys
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import pandas as pd
from config import DATA

sc = pd.read_csv(DATA / "s2_scene_scores.csv")
sc["date"] = pd.to_datetime(sc["date"])
sc["score"] = 0.85 * sc.clear_Pontchartrain + 0.15 * sc.clear_Maurepas
sc["month"] = sc.date.dt.to_period("M")

ok = sc[sc.clear_Pontchartrain >= 0.5]
monthly = ok.sort_values("score").groupby("month").tail(1)
ev = ok[((ok.date >= "2019-01-01") & (ok.date <= "2019-09-30")) |
        ((ok.date >= "2021-07-01") & (ok.date <= "2021-10-31"))]

ins = pd.read_csv(DATA / "insitu_wqp_lake_stations.csv")
sdates = pd.to_datetime(ins.ActivityStartDate.unique())
mu_ok = sc[sc.clear_Pontchartrain >= 0.3]
near = mu_ok[mu_ok.date.apply(lambda d: (abs((sdates - d).days) <= 1).any())]
cont = mu_ok[(mu_ok.date >= "2020-12-01") & (mu_ok.date <= "2022-10-11") & (mu_ok.clear_Pontchartrain >= 0.5)]

sel = pd.concat([monthly.assign(reason="monthly"), ev.assign(reason="event"),
                 near.assign(reason="matchup"), cont.assign(reason="matchup")])
sel = (sel.groupby("id").agg(date=("date", "first"), rel_orbit=("rel_orbit", "first"),
                             baseline=("baseline", "first"), score=("score", "first"),
                             clear_Pontchartrain=("clear_Pontchartrain", "first"),
                             clear_Maurepas=("clear_Maurepas", "first"),
                             reason=("reason", lambda r: ";".join(sorted(set(r)))))
       .reset_index().sort_values("date"))
sel["monthly"] = sel.reason.str.contains("monthly")
sel.to_csv(DATA / "s2_selected_scenes.csv", index=False)
print("months with a monthly scene:", sel.monthly.sum(), "of 72")
print("total scenes:", len(sel))
print(sel.reason.value_counts())
