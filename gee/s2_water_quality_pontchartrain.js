/*
 * Lake Pontchartrain and Lake Maurepas: Sentinel-2 L2A water-quality proxies
 * Google Earth Engine (JavaScript Code Editor) version of the Python workflow.
 *
 * NOT RUN HERE; FOR THE AUTHOR TO RUN IN THE GEE CODE EDITOR.
 * None of the numbers in the README come from this script. Results from GEE will
 * differ slightly (GEE uses the harmonized collection, its own resampling and a
 * monthly median composite rather than one clearest scene per month).
 *
 * Author: Anamika Dristi, Ph.D.
 * Proxies only: turbidity (Dogliotti et al. 2015 algorithm, FNU-like, uncalibrated),
 * NDCI (chlorophyll-a proxy) and B2/B3 (CDOM proxy). Not calibrated concentrations.
 */

var aoi = ee.Geometry.Rectangle([-90.65, 29.99, -89.77, 30.40]);
Map.centerObject(aoi, 10);

// COPERNICUS/S2_SR_HARMONIZED already removes the -1000 offset of baseline 04.00+
var s2 = ee.ImageCollection('COPERNICUS/S2_SR_HARMONIZED')
  .filterBounds(aoi)
  .filterDate('2019-01-01', '2025-01-01')
  .filter(ee.Filter.eq('MGRS_TILE', '15RYP'))
  .filter(ee.Filter.lt('CLOUDY_PIXEL_PERCENTAGE', 85));

// SCL classes masked: 0 no data, 1 saturated, 3 cloud shadow, 8/9 cloud, 10 cirrus, 11 snow
function prep(img) {
  var scl = img.select('SCL');
  var ok = scl.neq(0).and(scl.neq(1)).and(scl.neq(3)).and(scl.neq(8))
              .and(scl.neq(9)).and(scl.neq(10)).and(scl.neq(11));
  var r = img.select(['B2', 'B3', 'B4', 'B5', 'B8A', 'B11']).divide(10000);
  var mndwi = r.normalizedDifference(['B3', 'B11']).rename('MNDWI');
  var water = mndwi.gt(0);
  var red = r.select('B4'), nir = r.select('B8A');
  var tRed = red.multiply(228.1).divide(ee.Image(1).subtract(red.divide(0.1641)));
  var tNir = nir.multiply(3078.9).divide(ee.Image(1).subtract(nir.divide(0.2112)));
  var w = red.subtract(0.05).divide(0.02).clamp(0, 1);
  var turb = tRed.multiply(ee.Image(1).subtract(w)).add(tNir.multiply(w)).rename('TURB');
  var ndci = r.normalizedDifference(['B5', 'B4']).rename('NDCI');
  var cdom = r.select('B2').divide(r.select('B3')).rename('CDOM_B2B3');
  return turb.addBands(ndci).addBands(cdom).addBands(mndwi)
    .updateMask(ok.and(water))
    .copyProperties(img, ['system:time_start']);
}
var wq = s2.map(prep);

// Static open-water mask: MNDWI > 0 in >= 90% of clear observations, 120 m inland buffer
var freq = s2.map(function (img) {
  var scl = img.select('SCL');
  var ok = scl.neq(0).and(scl.neq(3)).and(scl.neq(8)).and(scl.neq(9)).and(scl.neq(10));
  return img.normalizedDifference(['B3', 'B11']).gt(0).updateMask(ok);
}).mean();
var openWater = freq.gte(0.9).focalMin({radius: 120, units: 'meters'}).selfMask();

// Monthly median composites, 2019-2024
var months = ee.List.sequence(0, 71);
var monthly = ee.ImageCollection(months.map(function (m) {
  var start = ee.Date('2019-01-01').advance(m, 'month');
  return wq.filterDate(start, start.advance(1, 'month')).median()
    .updateMask(openWater).set('system:time_start', start.millis());
}));

// Lake-wide median time series (whole AOI open water; draw lake polygons to split lakes)
var chart = ui.Chart.image.series({
  imageCollection: monthly.select('TURB'), region: aoi, reducer: ee.Reducer.median(),
  scale: 40
}).setOptions({title: 'Monthly median turbidity proxy (FNU-like), open water', pointSize: 3});
print(chart);

// Bonnet Carre Spillway 2019 (openings 27 Feb - 11 Apr, 10 May - 27 Jul) and Ida (29 Aug 2021)
var vis = {min: 2, max: 100, palette: ['ffffe5', 'fee391', 'fe9929', 'cc4c02', '662506']};
var pre19 = wq.filterDate('2019-01-01', '2019-02-27').median().select('TURB');
var open19 = wq.filterDate('2019-05-10', '2019-07-28').median().select('TURB');
var preIda = wq.filterDate('2021-07-01', '2021-08-29').median().select('TURB');
var postIda = wq.filterDate('2021-08-30', '2021-09-13').median().select('TURB');
Map.addLayer(pre19.updateMask(openWater), vis, 'Turbidity before BCS 2019');
Map.addLayer(open19.updateMask(openWater), vis, 'Turbidity BCS opening 2, 2019');
Map.addLayer(preIda.updateMask(openWater), vis, 'Turbidity before Ida');
Map.addLayer(postIda.updateMask(openWater), vis, 'Turbidity 1-14 days after Ida');
Map.addLayer(openWater, {palette: ['0000ff']}, 'Open-water mask', false);

// Export a reference median turbidity map (2021-2024) to Google Drive
var ref = wq.filterDate('2021-01-01', '2025-01-01')
  .filter(ee.Filter.or(ee.Filter.date('2021-01-01', '2021-08-29'),
                       ee.Filter.date('2021-11-01', '2025-01-01')))
  .select('TURB').median().updateMask(openWater);
Export.image.toDrive({image: ref, description: 'turbidity_proxy_median_2021_2024',
  region: aoi, scale: 40, crs: 'EPSG:32615', maxPixels: 1e9});
