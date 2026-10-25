# Metadata

FGDC Content Standard for Digital Geospatial Metadata (CSDGM, FGDC-STD-001-1998) records for the main geospatial products of this repository.

| Record | Data file | Product |
|---|---|---|
| [`turbidity_proxy_median_2021_2024_40m_fgdc.xml`](turbidity_proxy_median_2021_2024_40m_fgdc.xml) | `outputs/rasters/turbidity_proxy_median_2021_2024_40m.tif` | Median turbidity proxy of 38 reference monthly scenes, 2021 to 2024 (GeoTIFF, 40 m) |
| [`ndci_median_2021_2024_40m_fgdc.xml`](ndci_median_2021_2024_40m_fgdc.xml) | `outputs/rasters/ndci_median_2021_2024_40m.tif` | Median NDCI chlorophyll-a proxy of the same 38 scenes (GeoTIFF, 40 m) |
| [`turbidity_proxy_20190513_bonnet_carre_opening2_40m_fgdc.xml`](turbidity_proxy_20190513_bonnet_carre_opening2_40m_fgdc.xml) | `outputs/rasters/turbidity_proxy_20190513_bonnet_carre_opening2_40m.tif` | Turbidity proxy on 13 May 2019, second Bonnet Carre Spillway opening (GeoTIFF, 40 m) |

## What each record contains

Each XML file has the full CSDGM sections: identification (citation, abstract, purpose, time period, status, bounding coordinates, theme and place keywords with thesauri, access and use constraints), data quality (attribute accuracy, logical consistency, completeness, positional accuracy, lineage with source citations and processing steps), spatial data organization, spatial reference, entity and attribute definitions, distribution and metadata reference. Extents, coordinate systems, cell sizes, value ranges, class counts and field lists were read from the data files themselves; methods and accuracy figures come from the scripts and the main README.

## Validation

* Every file is well formed XML.
* Every file validates against the official FGDC CSDGM XML Schema (`fgdc-std-001-1998.xsd` and its section modules, downloaded from https://www.fgdc.gov/schemas/metadata/).
* Every file passes the USGS metadata parser `mp` version 2.9.50 (built from the USGS source at https://geology.usgs.gov/tools/metadata/) with no errors.

The standard requires a mailing address and a voice telephone for each contact. Only the city and state are given; the postal code and telephone fields read "Not provided". Please use email for contact.

## Viewing

Open the XML in a browser or a metadata editor, or turn it into text or HTML with `mp`, for example `mp turbidity_proxy_median_2021_2024_40m_fgdc.xml -t record.txt`.

Contact: Anamika Dristi, anamikadristi@gmail.com
