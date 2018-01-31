//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Preliminary: Evaluating the Impact of Various Methods of Atmospheric Correction on Detecting Surface Water in the Sahel of West Africa
//Kelsey E. Herndon, Graduate Research Assistant, SERVIR Science Coordination Office, The University of Alabama in Huntsville
//Date created: 9/15/2017
//Last modified: 11/1/2017  
//Purpose: this script evaluates the impact of different atmospheric correction methods on detecting surface water in the Ferlo Region of Senegal. Landsat TOA 
//reflectance, Landsat USGS surface reflectance, and the COST surface reflectance are assessed. 
//Notes: COST function by K. Markert SERVIR SCO
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//define imports
var lstoaCollection = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT_TOA"),
    geometry = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-16.4794921875, 16.5677576857484],
          [-16.5069580078125, 15.892659274817095],
          [-15.53192138671875, 15.91114960943624],
          [-15.63079833984375, 16.56249250837488]]]),
    lssr = ee.Image("users/keh0023/LC082050492016121201T1_SR_COMP_2"),
    geometry2 = /* color: #71d679 */ee.Geometry.Polygon(
        [[[-16.56463623046875, 16.227860792047952],
          [-16.1883544921875, 16.238409119210182],
          [-16.171875, 16.449256456410083],
          [-16.53167724609375, 16.45715879614139]]]);

//define point to filter landsat toa scenes
var pt = /* color: #d63000 */ee.Geometry.Point([-15.785980224609375, 16.337272136282255]);

//filter landsat scenes
var lstoa = ee.Image(lstoaCollection
  .filterBounds(pt)
  .filterDate('2016-12-11', '2016-12-13')
  .first());

//display resulting Landsat ToA scene
var trueColor1 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3};
Map.addLayer(lstoa, trueColor1,'Landsat TOA true color');

//display true color Landsat surface reflectance scene
var trueColor2 = {bands: ['b4', 'b3', 'b2'], min: 0, max: 3000};
Map.addLayer(lssr, trueColor2,'Landsat surface reflectance true color');

/////////////////////////////////////////////////////////////////////////////
//define function for atmospheric correction. COST function modified from script created by Kel Markert at SERVIR SCO - kel.markert@nasa.gov
function costz(img) {
  // implements the methods from http://info.asprs.org/publications/pers/96journal/september/1996_sep_1025-1036.pdf
  
  var Esun = {'OLI':[2004.57,1820.75,1549.49, 951.76,247.55,85.46],};

  var sunAngle = img.metadata('SUN_ELEVATION','sunAngle').subtract(ee.Image(90));
  var esunVal = Esun['OLI'];
  
  var bands = ['B2','B3','B4','B5','B6','B7'];
  var obsDate = ee.Date(img.get('system:time_start'));

  var jan01 = ee.Date.fromYMD(obsDate.get('year'),1,1);
  var doy_index = obsDate.difference(jan01,'day').toUint16();
  
  var dsun2 = ee.Image(1).subtract(ee.Image(0.01672).multiply(ee.Image(0.9856).multiply(ee.Image(doy_index).subtract(ee.Image(4)))).cos()).pow(ee.Image(2));
  
  var sr = ee.Image(0);
  var piImg = ee.Image(Math.PI);
  
  for (var i=0;i<bands.length;i++){
    var esunImg = ee.Image(esunVal[i]);

    var band = img.select(bands[i]);
    var denom = esunImg.multiply(sunAngle.multiply(piImg.divide(180)).cos());
    var numer = piImg.multiply(dsun2);
    var radiance = band.multiply(denom).divide(numer).rename('rad');
    
    var darkObject = radiance.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:img.geometry(),scale:150,bestEffort:true});
    var doImg = ee.Image(ee.Number(darkObject.get('rad')));

    var temp = ee.Image(0.01).multiply(esunImg).multiply(sunAngle.multiply(piImg.divide(180)).cos()).divide(dsun2.multiply(piImg));
    var Lhaze = doImg.subtract(temp);

    var band_sr = piImg.multiply(radiance.subtract(Lhaze)).multiply(dsun2).divide(esunImg.multiply(sunAngle.multiply(piImg.divide(180)).cos()));
    
    sr = sr.addBands(band_sr.rename(bands[i]));
  }
  
  return sr.set('system:time_start',obsDate);
}
//////////////////////////////////////////////////////////////////////////

//apply COST function to Landsat TOA scene
var lsCOST = costz(lstoa);

/////////////////////////////////////////////////////////////////////////
//display results 
var trueColor1 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3};
var trueColor2 = {bands: ['b4', 'b3', 'b2'], min: 0, max: 3000};
//display COST atmospherically corrected scene (true color)
Map.addLayer(lsCOST, trueColor1, 'Landsat with COST atmospheric correction');

//////////////////////////////////////////////////////////////////////////
//calculate mndwi and threshold for Landsat COST 
var mndwiCOST = lsCOST.normalizedDifference(['B3', 'B6']);
var h20maskCOST = mndwiCOST.gt(-0.2);
//display mndwi 
var waterPalette = ['yellow', 'blue'];
Map.addLayer(mndwiCOST, {min: -1, max: 1, palette: waterPalette}, 'mndwi COST');
//display water mask
var waterMaskPalette = ['white','blue'];
Map.addLayer(h20maskCOST, {min: 0, max: 1, palette: waterMaskPalette}, 'COST Water Mask');

//calculate mndwi and threshold for Landsat surface reflectance 
var mndwisr = lssr.normalizedDifference(['b3', 'b6']);
var h20masksr = mndwisr.gt(-0.2);
//display mndwi in user defined palette
Map.addLayer(mndwisr, {min: -1, max: 1, palette: waterPalette}, 'mndwi surface reflectance');
//display water mask in user defined palette
var waterMaskPalette = ['white','blue'];
Map.addLayer(h20masksr, {min: 0, max: 1, palette: waterMaskPalette}, 'surface reflectance Water Mask');

//calculate mndwi and threshold for Landsat TOA reflectance
var mndwiTOA = lstoa.normalizedDifference(['B3', 'B6']);
var h2omaskTOA = mndwiTOA.gt(-0.2);
//display mndwi in user defined palette
Map.addLayer(mndwiTOA, {min: -1, max: 1, palette: waterPalette}, 'mndwi TOA');
//display water mask in user defined palette
var waterMaskPalette = ['white','blue'];
Map.addLayer(h2omaskTOA, {min: 0, max: 1, palette: waterMaskPalette}, 'TOA Water Mask');

//calculate difference between mndwi for Landsat COST and Landsat sr
var mndwiDiff = mndwisr.subtract(mndwiCOST);
var mndwiDiffPalette = ['red','orange', 'white', 'orange','red'];
Map.addLayer(mndwiDiff,{palette: mndwiDiffPalette}, 'mNDWI Difference (COST and SR'); 

//histogram including all three water masks
var options = {
  format: 'short',
  hAxis: {title: 'Water (0-1 is No Water, 1-2 is Water Present)'},
  vAxis: {title: 'Frequency'},
  title: 'Histogram of Water Classification'
};
var images = [h20maskCOST, h20masksr, h2omaskTOA];
var histogram2 = ui.Chart.image.histogram(
  images, geometry2).
  setSeriesNames(['COST Water Mask', 'ESPA Water Mask', 'TOA Water Mask'])
  .setOptions(options);
print (histogram2);
