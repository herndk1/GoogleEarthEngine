//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Preliminary: Spectral Unmixing Analysis in GEE to Assess Changes in water quality Resulting from Illegal Gold Mining in the Peruvian Amazon
//Kelsey E. Herndon, Graduate Research Assistant, SERVIR Science Coordination Office, The University of Alabama in Huntsville
//Date Created: 10/15/2017
//Last Modified: 11/3/2017 
//Purpose: The purpose of this script is to identify a pre- and post-mining scene in the Peruvian Amazon, define endmembers within the Landsat scene, apply
//a spectral unmixing analysis (linear), and measure the change in the endmember proportions at several "virtual gage" locations
//Notes: COST atmospheric correction function provided by K. Markert at SERVIR Science Coordination Office 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//define imports
var ls5 = ee.ImageCollection("LANDSAT/LT5_L1T_TOA"),
    ls8 = ee.ImageCollection("LANDSAT/LC8_L1T_TOA"),
    pt = /* color: #d63000 */ee.Geometry.Point([-70.3729248046875, -12.6927606739619]),
    waterL5AC = /* color: #bf04c2 */ee.Geometry.Polygon(
        [[[-70.40107727050781, -13.151200050788454],
          [-70.40223598480225, -13.151994055802723],
          [-70.40215015411377, -13.153791005025262],
          [-70.40120601654053, -13.155504362948442],
          [-70.39970397949219, -13.156716243029326],
          [-70.3998327255249, -13.15583867529331],
          [-70.40107727050781, -13.154208898305768],
          [-70.40124893188477, -13.152453741741114]]]),
    soilL5AC = /* color: #ff0000 */ee.Geometry.Polygon(
        [[[-70.40249347686768, -13.134232803135985],
          [-70.40605545043945, -13.137994116868239],
          [-70.40528297424316, -13.138244869068094],
          [-70.40369510650635, -13.136782144291482],
          [-70.40202140808105, -13.135152240702864],
          [-70.4018497467041, -13.133647704709983]]]),
    vegL5AC = /* color: #00ff00 */ee.Geometry.Polygon(
        [[[-70.54269790649414, -13.168583993295154],
          [-70.54149627685547, -13.174434081120769],
          [-70.52553176879883, -13.17727550190674],
          [-70.52450180053711, -13.170422607388248]]]),
    gage4 = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[-69.89055633544922, -13.038090833182336],
          [-69.88924741744995, -13.038864298464551],
          [-69.88894701004028, -13.038508922824336],
          [-69.89104986190796, -13.037547315591555]]]),
    gage1 = /* color: #bf04c2 */ee.Geometry.Polygon(
        [[[-69.97280359268188, -13.01969418982366],
          [-69.97454166412354, -13.01971509587577],
          [-69.97469186782837, -13.020154122562499],
          [-69.97417688369751, -13.019945062332502]]]),
    gage2 = /* color: #ff0000 */ee.Geometry.Polygon(
        [[[-69.92169141769409, -13.036481181465001],
          [-69.92351531982422, -13.036543895364307],
          [-69.9234938621521, -13.036920178426369],
          [-69.92154121398926, -13.03689927382682]]]),
    gage3 = /* color: #00ff00 */ee.Geometry.Polygon(
        [[[-69.90137100219727, -13.033679944420493],
          [-69.90173578262329, -13.033679944420493],
          [-69.90171432495117, -13.034098041514747],
          [-69.90119934082031, -13.034265280154692],
          [-69.90074872970581, -13.034578852299996]]]),
    gage5 = /* color: #0000ff */ee.Geometry.Polygon(
        [[[-69.90222930908203, -13.041979039583738],
          [-69.90291595458984, -13.042188081207465],
          [-69.90321636199951, -13.042689780383581],
          [-69.90261554718018, -13.042522547437951]]]),
    gage6 = /* color: #999900 */ee.Geometry.Polygon(
        [[[-69.8682188987732, -13.042020847922622],
          [-69.868905544281, -13.042438930922719],
          [-69.86817598342896, -13.042334410238956]]]),
    gage7 = /* color: #009999 */ee.Geometry.Polygon(
        [[[-69.85972166061401, -13.047351353232578],
          [-69.86042976379395, -13.047330449514606],
          [-69.86032247543335, -13.047644005098542],
          [-69.85982894897461, -13.047664908789994]]]);

////////////////////////////////////////////////////////
//select post image
var post = ee.Image(ls8
   .filterBounds(pt)
   .filterDate('2013-01-01', '2016-09-27')
   .sort('CLOUD_COVER')
   .first());

//display POST Landsat 8 TOA reflectance
var trueColor1 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.5};
Map.addLayer(post,trueColor1,'Post Original');

Map.centerObject(post,11);

/////////////////////////////////////////////////////////////////////
//apply COST atmospheric correction to post image: 
//define function for atmospheric correction 
function costz8(img) {
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
var postAC = costz8(post);

//add post, true color, atmospherically corrected image to map
Map.addLayer(postAC, trueColor1, 'Post Image with Atmospheric Correction');

/////////////////////////////////////////////////////////////////////////
//process pre image
//select pre image
var pre = ee.Image(ls5
   .filterBounds(pt)
   .filterDate('2005-01-01', '2010-12-31')
   .sort('CLOUD_COVER')
   .first());
var trueColor2 =  {bands: ['B3', 'B2', 'B1'], min: 0, max: 0.5};
Map.addLayer(pre,trueColor2, 'Pre Original');

/////////////////////////////////////////////////////////////////////
//apply atmospheric correction to pre image 
//define function for atmospheric correction 
function costz5(img) {
  // implements the methods from http://info.asprs.org/publications/pers/96journal/september/1996_sep_1025-1036.pdf
  
  var Esun = {'TM':[1958, 1827, 1551, 1036, 214.9, 80.65],};//Esun values taken from: https://landsat.usgs.gov/esun
  
  var sunAngle = img.metadata('SUN_ELEVATION','sunAngle').subtract(ee.Image(90));
  var esunVal = Esun['TM'];
  
  var bands = ['B1', 'B2','B3','B4','B5','B7'];
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

//apply atmospheric correction to pre image 
var preAC = costz5(pre);

//display pre true color, atmospherically corrected image

Map.addLayer(preAC, trueColor2, 'Pre Image with Atmospheric Correction');

//////////////////////////////////////////////////////////////////////////////////////////////
//Get values of spectral endmemebers - these are taken from pre-image (Landsat 5) after atmospheric correction image

var soilDict = preAC.reduceRegion(ee.Reducer.mean(),soilL5AC, 30);
var soilList = soilDict.values(['B1', 'B2', 'B3','B4', 'B5', 'B7']);
print (soilDict);
print (soilList);

var vegDict = preAC.reduceRegion(ee.Reducer.mean(),vegL5AC, 30);
var vegList = vegDict.values(['B1', 'B2', 'B3','B4', 'B5', 'B7']);
print (vegDict);
print (vegList);

var waterDict = preAC.reduceRegion(ee.Reducer.mean(),waterL5AC, 30);
var waterList = waterDict.values(['B1', 'B2', 'B3','B4', 'B5', 'B7']);
print (waterDict);
print (waterList);

///////////////////////////////////////////////////////////////////////
//select the post bands to be used in spectral unmixing
var postACb = postAC.select('B2', 'B3','B4', 'B5', 'B6','B7');
// Unmix post image.
var fractionsPost = postACb.unmix([soilList, vegList, waterList], true,true);
Map.addLayer(fractionsPost, {min:0, max: 0.6}, 'Post Unmixed Image');

//select pre bands for spectral unmixing
var preACb = preAC.select('B1', 'B2', 'B3','B4', 'B5', 'B7');
//unmix pre image
var fractionsPre = preACb.unmix([soilList, vegList, waterList],true,true);
Map.addLayer(fractionsPre, {min:0, max: 0.6}, 'Pre Unmixed Image');

////////////////////////////////////////////////////////////////////////////
//calculate ratio of "soil" to "water" for pre and post image
var preRatio = fractionsPre.select('band_0').divide(fractionsPre.select('band_2'));
Map.addLayer(preRatio,{min:-1, max:2},'PreRatio');

var postRatio = fractionsPost.select('band_0').divide(fractionsPost.select('band_2'));
Map.addLayer(postRatio, {min:-1, max: 2}, 'PostRatio');

//calculate difference in soil:water ratios
var diffRatio = postRatio.subtract(preRatio);
Map.addLayer(diffRatio, {min:-1, max:1}, 'DiffRatio');

//create an image where areas with an increase in the ratio of soil to water are given a value of 1 and 
//...areas where there was a decrease in the ratio of soil to water are given a value of 0
var diffRatioMask = diffRatio.gt(0);
Map.addLayer(diffRatioMask,{}, 'Mask');

///////////////////////////////////////////////////////////////////////////////////////////////////////////////
//identify points where water is present in both pre and post images

//calculate the Modified Normalized Difference Vegetation Index (MNDWI)using Landsat 8 Band 3 and Band 6
var mndwiPre = preAC.normalizedDifference(['B2', 'B5']);
var mndwiPost = postAC.normalizedDifference(['B3','B6']);

//create a water mask where MNDWI values greater than a threshold (in this case -0.1) are given a value 
//of 1 and values less than the threshold are given a value of zero
var h20maskPre = mndwiPre.gt(0.0);
var h20maskPost = mndwiPost.gt(0.0);

//display mndwi in user defined palette
var waterPalette = ['yellow', 'blue'];
Map.addLayer(mndwiPre, {min: -1, max: 1, palette: waterPalette}, 'Pre MNDWI');
Map.addLayer(mndwiPost, {min: -1, max: 1, palette: waterPalette}, 'Post MNDWI');

//display water mask in user defined palette
var waterMaskPalette = ['white','blue'];
Map.addLayer(h20maskPre, {min: 0, max: 1, palette: waterMaskPalette}, 'Pre Water Mask');
Map.addLayer(h20maskPost, {min: 0, max: 1, palette: waterMaskPalette}, 'Post Water Mask');

//create layer where 2 indicates water is present in both years
var permWater = h20maskPre.add(h20maskPost)
Map.addLayer(permWater,{},'Permanant Water')

//////////////////////////////////////////////////////////

///////////////////////////////////////////////////////////////////////////////
//display "virtual gage" values
var geometry2 = [gage1, gage2, gage3, gage4, gage5, gage6, gage7];
var images = [preRatio, postRatio];
var images2 = [fractionsPre.select('band_0'), fractionsPost.select('band_0')]
var options = {
  format: 'bar',
  hAxis: {title: 'Gage Number'}, 
  vAxis: {title: 'Ratio of Soil to Water'},
  title: 'Gage Measurements'
};

var chart2 = ui.Chart.image.byRegion(
  images2, geometry2,ee.Reducer.mean()).
  setSeriesNames(['Pre Soil Signature', 'Post Soil Signature'])
  .setOptions(options);
print (chart2);
