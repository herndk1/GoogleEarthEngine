///////////////////////////////////////////////////////////////////////////////////////////////
//Kelsey E. Herndon, Graduate Research Assistant, SERVIR Science Coordination Office, The 
//University of Alabama in Huntsville, kelsey.e.herndon@nasa.gov
//Date Created: 1/10/2017
//Last Update: 1/10/2017 
//Purpose: The purpose of this script is to apply canny edge detection and Otsu thresholding to multiple indices to compare
//their ability to identify surface water bodies in Niger
//Notes: This script is based on work by Nick Clinton (Google) and Gennadii Donchyts, http://www.mdpi.com/2072-4292/8/5/386
//The changes that were made include: limiting analysis to one scene,changing study area to West Africa, introducing 
//an additional for-loop that will apply the algorithm to a list of indices,
//new notes and parameters,and an additional raster that calculates the number of methods where a pixel is identified as water
///////////////////////////////////////////////////////////////////////////////////////////////

//define imports
var ls8 = ee.ImageCollection("LANDSAT/LC8_L1T_TOA"),
    pt = /* color: #d63000 */ee.Geometry.Point([5.1656341552734375, 15.389473691090894]);
//define image
var img = ee.Image(ls8
  .filterBounds(pt)
  .filterDate('2013-12-11', '2017-01-10')
  .sort('CLOUD_COVER')
  .first());
  
//display true color landsat image 
var TrueColorLS8 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.5};
Map.addLayer(img, TrueColorLS8, 'True Color Landsat Image')

//define the otsu function:///////////////////////////////////////////////////////////////// 
//otsu function is used in the detect water function
// Return the value that maximizes interclass variance (in the region).
var otsu = function(histogram) {
  var counts = ee.Array(ee.Dictionary(histogram).get('histogram'));
  var means = ee.Array(ee.Dictionary(histogram).get('bucketMeans'));
  var size = means.length().get([0]);
  var total = counts.reduce(ee.Reducer.sum(),[0]).get([0]);
  var sum = means.multiply(counts).reduce(ee.Reducer.sum(), [0]).get([0]);
  var mean = sum.divide(total);
  
  var indices = ee.List.sequence(1, size);
  //print (indices)
  
// Compute between sum of squares (BSS), where each mean partitions the data.
  var bss = indices.map(function(i) {
    var aCounts = counts.slice(0, 0, i);
    var aCount = aCounts.reduce(ee.Reducer.sum(), [0]).get([0]);
    var aMeans = means.slice(0, 0, i);
    var aMean = aMeans.multiply(aCounts)
        .reduce(ee.Reducer.sum(), [0]).get([0])
        .divide(aCount);
    var bCount = total.subtract(aCount);
    var bMean = sum.subtract(aCount.multiply(aMean)).divide(bCount);
    return aCount.multiply(aMean.subtract(mean).pow(2)).add(
           bCount.multiply(bMean.subtract(mean).pow(2)));
  });
  
  print(ui.Chart.array.values(ee.Array(bss), 0, means));
  
  // Return the mean value corresponding to the maximum BSS.
  return means.sort(bss).get([-1]);
};

/////////////////////////////////////////////////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////////////////////////////////////////////////

//define detect water function
function detectWater(image, colour, name) {

    // *compute MNDWI with SWIR 1
  var mndwi1 = image.normalizedDifference(MNDWI1_BANDS).set({'NAME':'mndwi1'})
  var mndwi1_min = -0.3
  var mndwi1_max = 0.6
  var mndwi1_vis = {min: mndwi1_min, max: mndwi1_max}
  Map.addLayer(mndwi1, mndwi1_vis, 'MNDWI1 (B/W)', false)
  
  //*compute MNDWI with SWIR 2
  var mndwi2 = image.normalizedDifference(MNDWI2_BANDS).set({'NAME':'mndwi2'})
  var mndwi2_min = -0.3
  var mndwi2_max = 0.6
  var mndwi2_vis = {min: mndwi2_min, max: mndwi2_max}
  Map.addLayer(mndwi2, mndwi2_vis, 'MNDWI2 (B/W)', false)
  
  //*compute AWEISH 
  var aweiSH = image.expression('BLUE + 2.5 * GREEN - 1.5 * (NIR + SWIR1) - 0.25 * SWIR2', {
  'BLUE': img2.select('blue'),
  'GREEN': img2.select('green'),
  'NIR': img2.select('nir'),
  'SWIR1': img2.select('swir1'),
  'SWIR2': img2.select('swir2')
  }).rename('nd');
  var aweiSH_min = -1
  var aweiSH_max = 1
  var aweiSH_vis = {min: aweiSH_min, max: aweiSH_max}
  Map.addLayer(aweiSH, aweiSH_vis, 'AWEI (B/W)', false)
  
  //*create a list of the rasters created from the above indices
  var waterIndices = [mndwi1, mndwi2, aweiSH]
  var names = ['mndwi1','mndwi2', 'aweiSH']
  
  //*apply canny edge detection and otsu thresholding for all indices in waterIndices list
  for (var i=0;i<waterIndices.length; i++){
  
    var index = waterIndices[i];
    var NAME = names[i];
    print (NAME)

    
    // detect sharp changes in NDWI
    var canny = ee.Algorithms.CannyEdgeDetector(index.clip(bounds), 0.7, 0.5);
    canny = canny.mask(canny).clip(bounds)
  
    Map.addLayer(canny, {min: 0, max: 1, palette: 'FF0000'}, 'canny '+ NAME , false);
  
    // buffer around NDWI edges
    var cannyBuffer = canny.focal_max(ee.Number(scale).multiply(1.5), 'square', 'meters');
      var index_canny = index.mask(cannyBuffer)
    
    Map.addLayer(cannyBuffer, {}, 'canny Buffer ' + NAME)
  
    if(Map.getScale() > 400) {
      throw 'Error: zoom in to at least 1km scale to apply dynamic thresholding'
    }

    // print charts
    print(Chart.image.histogram(index, ee.Geometry(bounds), scale, 255).setOptions({title: NAME, vAxis: { gridlines: { count: 0 } }, hAxis: { gridlines: { count: 0 }, viewWindow:{max:mndwi1_min, min:mndwi1_max} }}));
    print(Chart.image.histogram(index_canny, ee.Geometry(bounds), scale, 255).setOptions({title: NAME + ' around canny', vAxis: { gridlines: { count: 0 } }, hAxis: { gridlines: { count: 0 },viewWindow:{max:mndwi1_min, min:mndwi1_max} }}));
    Map.addLayer(index_canny, mndwi1_vis, NAME + ' around canny', false);
  
    // simple 0 thresholding
    var VIS_WATER = {min:0.03, max:colour} //max is [0.4,0.4,0.5]
    var water0 = index.gt(0)
    Map.addLayer(image.mask(water0), VIS_WATER, 'water (0)' + NAME, false)
  
    // compute threshold using Otsu thresholding
    var hist = index_canny.reduceRegion(ee.Reducer.histogram(255), bounds, scale)
    var index_threshold = otsu(hist.get('nd'));//'nd'
    index_threshold = ee.Number(mndwi1_max).min(ee.Number(mndwi1_min).max(index_threshold))
    print(NAME, 'Detected threshold' , index_threshold)

    // show water mask
    var water = index.gt(index_threshold)
    Map.addLayer(image.mask(water), VIS_WATER, 'water ' + NAME);
  
    // show edge around water mask
    var canny = ee.Algorithms.CannyEdgeDetector(water, 0.99, 0.3);
    canny = canny.mask(canny)
    Map.addLayer(canny, {palette:'aaaaff'}, 'water (boundary)', false)
    Map.addLayer(canny, {palette:'000000'}, 'water (boundary) black', false)
  }
}

//////////////////////////////////////////////////////////////////////////////////////////////////

// Parameters
var bounds = ee.Geometry(Map.getBounds(true)).buffer(-100);

var VIS_IMAGE = {min:0.03, max:[0.4,0.4,0.5]}
var MNDWI1_BANDS = ['green', 'swir1']
var MNDWI2_BANDS =['green', 'swir2']
var scale = 10

// Bands
var l8_bands = ['B2','B3','B4','B5','B6','B7']
var bandNames = ['blue','green','red','nir','swir1','swir2']

var img2 = img.select(l8_bands, bandNames)

////////////////////////////////////////////////////////////////////////////////////
//Apply function to image
detectWater(img2, 0.1, 2015)
