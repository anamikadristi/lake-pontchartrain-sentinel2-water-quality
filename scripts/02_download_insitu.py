"""Step 2: download public in-situ water-quality data for the two lakes.

(a) USGS continuous turbidity (63680, FNU) and chlorophyll fluorescence (32316, ug/L)
    at USGS 300138090064700, New Canal Lighthouse at New Orleans, LA, from the keyless
    USGS Water Data OGC API (api.waterdata.usgs.gov). The legacy NWIS site service
    was used only to find which stations in the area report turbidity.
(b) Discrete turbidity and chlorophyll-a results from the Water Quality Portal
    (waterqualitydata.us), mainly LDEQ ambient monitoring, 2019-2024.
Raw files go to WORK/insitu; small cleaned tables go to data/.
"""
import sys, io, time, json
from pathlib import Path
sys.path.insert(0, str(Path(__file__).resolve().parent))
import requests
import pandas as pd
from config import WORK, DATA, START, END

RAW = WORK / "insitu"
RAW.mkdir(parents=True, exist_ok=True)
SITE = "USGS-300138090064700"
OGC = "https://api.waterdata.usgs.gov/ogcapi/v0/collections/continuous/items"


def get_continuous(pcode):
    out = RAW / f"usgs_continuous_{pcode}.csv"
    if out.exists():
        return pd.read_csv(out)
    recs = []
    # the API limits one query to 1100 days, so request one calendar year at a time
    for yr in range(int(START[:4]), int(END[:4]) + 1):
        url = OGC
        params = dict(f="json", monitoring_location_id=SITE, parameter_code=pcode,
                      time=f"{yr}-01-01T00:00:00Z/{yr}-12-31T23:59:59Z", limit=10000,
                      properties="time,value,unit_of_measure,approval_status,qualifier,statistic_id")
        while url:
            for attempt in range(5):
                r = requests.get(url, params=params, timeout=180)
                if r.status_code == 200:
                    break
                time.sleep(20 * (attempt + 1))
            r.raise_for_status()
            js = r.json()
            recs += [f["properties"] for f in js["features"]]
            nxt = [l["href"] for l in js.get("links", []) if l.get("rel") == "next"]
            url, params = (nxt[0], None) if nxt and js["features"] else (None, None)
        print(pcode, yr, len(recs))
    df = pd.DataFrame(recs)
    df.to_csv(out, index=False)
    return df


for pc, name in [("63680", "turbidity_fnu"), ("32316", "chl_fluor_ugL")]:
    d = get_continuous(pc)
    d["time"] = pd.to_datetime(d["time"], utc=True)
    d["value"] = pd.to_numeric(d["value"], errors="coerce")
    print(name, len(d), d.time.min(), d.time.max())

# ---- Water Quality Portal discrete data
WQP = "https://www.waterqualitydata.us/data/Result/search"
wq_out = RAW / "wqp_results.csv"
if not wq_out.exists():
    params = {"bBox": "-90.65,29.99,-89.72,30.40",
              "characteristicName": ["Turbidity", "Chlorophyll a"],
              "startDateLo": "01-01-2019", "startDateHi": "12-31-2024",
              "mimeType": "csv", "zip": "no", "dataProfile": "resultPhysChem"}
    r = requests.get(WQP, params=params, timeout=600)
    r.raise_for_status()
    wq_out.write_bytes(r.content)
st_out = RAW / "wqp_stations.csv"
if not st_out.exists():
    params = {"bBox": "-90.65,29.99,-89.72,30.40",
              "characteristicName": ["Turbidity", "Chlorophyll a"],
              "startDateLo": "01-01-2019", "startDateHi": "12-31-2024", "mimeType": "csv"}
    r = requests.get("https://www.waterqualitydata.us/data/Station/search", params=params, timeout=600)
    r.raise_for_status()
    st_out.write_bytes(r.content)

wq = pd.read_csv(wq_out, low_memory=False)
st = pd.read_csv(st_out, low_memory=False)
print("WQP raw results:", len(wq), "stations:", len(st))
st = st[["MonitoringLocationIdentifier", "MonitoringLocationName", "MonitoringLocationTypeName",
         "LatitudeMeasure", "LongitudeMeasure"]]
# Keep open-lake stations only (names starting with Lake Pontchartrain / Lake Maurepas).
lake = st[st.MonitoringLocationName.str.match(r"^\s*Lake (Pontchartrain|Maurepas)", na=False)]
wq = wq.drop(columns=[c for c in lake.columns if c in wq.columns and c != "MonitoringLocationIdentifier"])
wq = wq.merge(lake, on="MonitoringLocationIdentifier", how="inner")
wq["value"] = pd.to_numeric(wq["ResultMeasureValue"], errors="coerce")
wq = wq.dropna(subset=["value"])
keep = ["OrganizationIdentifier", "MonitoringLocationIdentifier", "MonitoringLocationName",
        "LatitudeMeasure", "LongitudeMeasure", "ActivityStartDate", "ActivityStartTime/Time",
        "ActivityStartTime/TimeZoneCode", "ActivityDepthHeightMeasure/MeasureValue",
        "CharacteristicName", "value", "ResultMeasure/MeasureUnitCode",
        "ResultSampleFractionText", "ResultAnalyticalMethod/MethodName"]
keep = [k for k in keep if k in wq.columns]
wq = wq[keep].sort_values(["CharacteristicName", "ActivityStartDate"])
wq.to_csv(DATA / "insitu_wqp_lake_stations.csv", index=False)
print(wq.groupby(["CharacteristicName", "ResultMeasure/MeasureUnitCode"]).agg(
    n=("value", "size"), stations=("MonitoringLocationIdentifier", "nunique"),
    first=("ActivityStartDate", "min"), last=("ActivityStartDate", "max")))
