//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Preliminary: Spectral Unmixing Analysis in GEE to Assess the Changes water quality from Gold Mining in the Peruvian Amazon
//Kelsey E. Herndon, Graduate Research Assistant, SERVIR Science Coordination Office, The University of Alabama in Huntsville
//Date Created: 10/15/2017
//Last Modified: 11/3/2017
//Purpose: The purpose of this script is to identify a pre- and post-mining Landsat scene covering the Peruvian Amazaon, define endmembers within the Landsat scene, apply
//spectral unmixing analysis, and measure the relative change in the soil percentage of river pixels at several "virtual gage" locations
//Notes: COST atmospheric correction function modified from original by Kel Markert, SERVIR SCO; 
//BRDF correction function modified from original by Daniel Wiell and Erik Lindquist of the United Nations Food and Agriculture Organization
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//def imports
var pt = /* color: #d63000 */ee.Geometry.Point([-70.3729248046875, -12.6927606739619]),
    ls5 = ee.ImageCollection("LANDSAT/LT5_L1T"),
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
          [-69.97415542602539, -13.01997642137824],
          [-69.97448801984706, -13.019960741977744],
          [-69.97472405433655, -13.02007049849167],
          [-69.97440218928193, -13.020049592667087],
          [-69.97408032417297, -13.020112310530608]]]),
    gage2 = /* color: #ff0000 */ee.Geometry.Polygon(
        [[[-69.92156267166138, -13.03673203696687],
          [-69.92253899562724, -13.036773847785353],
          [-69.92350459098816, -13.03674248927392],
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
          [-69.90270137786865, -13.042386670586366]]]),
    gage6 = /* color: #999900 */ee.Geometry.Polygon(
        [[[-69.8682188987732, -13.042020847922622],
          [-69.868905544281, -13.042438930922719],
          [-69.86842274665833, -13.042365766448723]]]),
    gage7 = /* color: #009999 */ee.Geometry.Polygon(
        [[[-69.85972166061401, -13.047351353232578],
          [-69.86042976379395, -13.047330449514606],
          [-69.86032247543335, -13.047644005098542],
          [-69.85982894897461, -13.047664908789994]]]),
    soilL5ACl8 = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-70.40502548217773, -13.136907521042405],
          [-70.40614128112793, -13.137994116868239],
          [-70.40528297424316, -13.138244869068094]]]),
    ls8 = ee.ImageCollection("LANDSAT/LC8_L1T"),
    gage8 = /* color: #98ff00 */ee.Geometry.Polygon(
        [[[-69.97559309005737, -13.020321370619492],
          [-69.97492790222168, -13.020321370619492],
          [-69.97490644454956, -13.019986874392611],
          [-69.97555017471313, -13.020007780420013]]]),
    gage9 = /* color: #0b4a8b */ee.Geometry.Polygon(
        [[[-69.92850422859192, -13.037035153692301],
          [-69.92769956588745, -13.037076964424347],
          [-69.92694854736328, -13.037087415159164],
          [-69.92681980133057, -13.03694108302415],
          [-69.92860078811646, -13.03694108302415]]]),
    gage10 = /* color: #ffc82d */ee.Geometry.Polygon(
        [[[-69.89656448364258, -13.040139465677369],
          [-69.89720821380615, -13.040682977573535],
          [-69.8969292640686, -13.040703881853387]]]),
    gage11 = /* color: #00ffff */ee.Geometry.Polygon(
        [[[-69.87786412239075, -13.044989221919154],
          [-69.87802505493164, -13.045020577792341],
          [-69.87802505493164, -13.045490915412893]]]),
    gage12 = /* color: #bf04c2 */ee.Geometry.Polygon(
        [[[-69.8615026473999, -13.038916559545054],
          [-69.8617172241211, -13.038697062932737],
          [-69.86145973205566, -13.039219673594314]]]),
    gage13 = /* color: #ff0000 */ee.Geometry.Polygon(
        [[[-69.85755443572998, -13.047246834625069],
          [-69.85826253890991, -13.04748722735628],
          [-69.85850930213928, -13.047717168010902],
          [-69.8580050464036, -13.047581294586161]]]),
    gage14 = /* color: #00ff00 */ee.Geometry.Polygon(
        [[[-69.84025955200195, -13.048877319870359],
          [-69.84086036682129, -13.049379005474716],
          [-69.8406457901001, -13.050465987458004]]]),
    gage15 = /* color: #0000ff */ee.Geometry.Polygon(
        [[[-69.84738349914551, -13.051490254567963],
          [-69.84739422798157, -13.052159161409856],
          [-69.84690070152283, -13.053078905361414]]]),
    gage16 = /* color: #999900 */ee.Geometry.Polygon(
        [[[-69.840989112854, -13.043567751493024],
          [-69.84137535095215, -13.043860407836902],
          [-69.84109640195283, -13.044351650532963],
          [-69.84090328216553, -13.04540729990129]]]),
    gage17 = /* color: #009999 */ee.Geometry.Polygon(
        [[[-69.83385443687439, -13.042438930922719],
          [-69.83367204666138, -13.042135820818096],
          [-69.8352599143982, -13.041184679802411]]]),
    gage18 = /* color: #ff00ff */ee.Geometry.Polygon(
        [[[-69.8291015625, -13.045783569480854],
          [-69.82919812202454, -13.045302780472015],
          [-69.83026027679443, -13.044696566911108],
          [-69.82938051351084, -13.045396848007151]]]),
    gage19 = /* color: #ff9999 */ee.Geometry.Polygon(
        [[[-69.79820251464844, -13.02613317046672],
          [-69.79773044586182, -13.02613317046672],
          [-69.79901790618896, -13.025296948671416]]]),
    gage20 = /* color: #99ff99 */ee.Geometry.Polygon(
        [[[-69.77725982666016, -13.01802169994064],
          [-69.77833271026611, -13.019025195225034],
          [-69.77773189544678, -13.01889975853671]]]),
    gage21 = /* color: #9999ff */ee.Geometry.Polygon(
        [[[-70.00183582305908, -13.011749762361795],
          [-70.00247955322266, -13.01187520266847],
          [-70.0019645690918, -13.012167896470544]]]),
    gage22 = /* color: #ffff99 */ee.Geometry.Polygon(
        [[[-69.99930381774902, -13.012418776597308],
          [-69.9994969367981, -13.012209709842631],
          [-69.9997329711914, -13.012105176571923],
          [-69.9999475479126, -13.012000642911662],
          [-69.99969005611854, -13.012397870039935]]]),
    gage23 = /* color: #99ffff */ee.Geometry.Polygon(
        [[[-69.97982025146484, -13.016683699907208],
          [-69.97991681085637, -13.016913669231045],
          [-69.98008847236633, -13.017028653731918],
          [-69.97978806495667, -13.017049560008903]]]),
    gage25 = /* color: #d63000 */ee.Geometry.Polygon(
        [[[-69.9695634841919, -13.0142794629469],
          [-69.96974587440491, -13.014488527956138],
          [-69.96986389160156, -13.014781218672844],
          [-69.96962785705625, -13.014613966793284]]]);

var vizParams = {bands: 'red, green, blue', min: '200, 400, 600', max: '2400, 2200, 2400', gamma: 1.2}
var vizParamsFalse = {bands: 'nir, swir1, red', min: 0, max: 5000, gamma: 1.7}

// ************** define pre and post images **************
//define pre image and display
var pre = ee.Image(ls5
   .filterBounds(pt)
   .filterDate('2009-01-01', '2010-12-31')
   .sort('CLOUD_COVER')
   .first());
var trueColor2 =  {bands: ['B3', 'B2', 'B1'], min: 0, max: 150};
Map.addLayer(pre,trueColor2, 'Pre Original');


//define post image  and display
var post = ee.Image(ls8
   .filterBounds(pt)
   .filterDate('2013-01-01', '2016-09-27')
   .sort('CLOUD_COVER')
   .first());
var trueColor1 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 25000};
Map.addLayer(post,trueColor1,'Post Original');

//center map on Landsat scene
Map.centerObject(post,11);

//////////////////////////////////////////////////////////////////////////////////////////
//rename bands in pre and post scenes to common names 
var pre2 = pre.select(['B1','B2','B3','B4','B5','B7'],['blue','green','red','nir','swir1','swir2']);
var post2 = post.select(['B2','B3','B4','B5','B6','B7'],['blue','green','red','nir','swir1','swir2']);


//////////////////////////////////////////////////////////////////////////////////////////
//convert from DN to TOA Radiance for pre scene (Landsat 5)
var bandsls5 = ['B1','B2','B3','B4','B5','B7'];

//convert DN to Radiance
//define variables 
//make fuction
//
var B1Mls5 = ee.Image(pre.get('RADIANCE_MULT_BAND_1'));
var B1Als5 = ee.Image(pre.get('RADIANCE_ADD_BAND_1'));
print ('B1M', B1Mls5);
print ('B1A', B1Als5);

var B2Mls5 = ee.Image(pre.get('RADIANCE_MULT_BAND_2'));
var B2Als5 = ee.Image(pre.get('RADIANCE_ADD_BAND_2'));
print ('B2M', B2Mls5);
print ('B2A', B2Als5);

var B3Mls5 = ee.Image(pre.get('RADIANCE_MULT_BAND_3'));
var B3Als5 = ee.Image(pre.get('RADIANCE_ADD_BAND_3'));
print ('B3M', B3Mls5);
print ('B3A', B3Als5); 

var B4Mls5 = ee.Image(pre.get('RADIANCE_MULT_BAND_4'));
var B4Als5 = ee.Image(pre.get('RADIANCE_ADD_BAND_4'));
print ('B4M', B4Mls5);
print ('B4A', B4Als5);

var B5Mls5 = ee.Image(pre.get('RADIANCE_MULT_BAND_5'));
var B5Als5 = ee.Image(pre.get('RADIANCE_ADD_BAND_5'));
print ('B5M', B5Mls5);
print ('B5A', B5Als5);

var B7Mls5 = ee.Image(pre.get('RADIANCE_MULT_BAND_7'));
var B7Als5 = ee.Image(pre.get('RADIANCE_ADD_BAND_7'));
print ('B7M', B7Mls5);
print ('B7A', B7Als5);

var BMls5 = [B1Mls5, B2Mls5, B3Mls5, B4Mls5, B5Mls5, B7Mls5];
var BAls5 = [B1Als5, B2Als5, B3Als5, B4Als5, B5Als5, B7Als5];
print (BMls5);
print (BAls5);

//create an empty ee.image to save new radiance bands to
var ls5toaRad = ee.Image(0);

//calculate radiance by looping through bands
for (var i=0;i<bandsls5.length; i++){
  
  var band = pre.select(bandsls5[i]);
  print (band);
  var mult = ee.Number(BMls5[i]);
  print(mult);
  var add = ee.Number(BAls5[i]);
  print(add);
  var Radls5 = (band.multiply(ee.Image(mult)).add(ee.Image(add)));
  
  ls5toaRad = ls5toaRad.addBands(Radls5.rename(bandsls5[i]));
}  

//check result
print (ls5toaRad);

//display result
var trueColor3 = {bands: ['B3', 'B2', 'B1'], min: 0, max: 100};
Map.addLayer(ls5toaRad, trueColor3,'Pre Landsat Radiance True Color');  

/////////////////////////////////////////////////////////////////////////////////////////
//convert from DN to TOA Radiance for post scene (Landsat 8) 
var bandsls8 = ['B2','B3','B4','B5','B6','B7'];

//convert DN to Radiance
//define variables 

var B2Mls8 = ee.Image(post2.get('RADIANCE_MULT_BAND_2'));
var B2Als8 = ee.Image(post2.get('RADIANCE_ADD_BAND_2'));
print ('B2M', B2Mls8);
print ('B2A', B2Als8);

var B3Mls8 = ee.Image(post2.get('RADIANCE_MULT_BAND_3'));
var B3Als8 = ee.Image(post2.get('RADIANCE_ADD_BAND_3'));
print ('B3M', B3Mls8);
print ('B3A', B3Als8);

var B4Mls8 = ee.Image(post2.get('RADIANCE_MULT_BAND_4'));
var B4Als8 = ee.Image(post2.get('RADIANCE_ADD_BAND_4'));
print ('B4M', B4Mls8);
print ('B4A', B4Als8);

var B5Mls8 = ee.Image(post2.get('RADIANCE_MULT_BAND_5'));
var B5Als8 = ee.Image(post2.get('RADIANCE_ADD_BAND_5'));
print ('B5M', B5Mls8);
print ('B5A', B5Als8);

var B6Mls8 = ee.Image(post2.get('RADIANCE_MULT_BAND_6'));
var B6Als8 = ee.Image(post2.get('RADIANCE_ADD_BAND_6'));
print ('B6M', B6Mls8);
print ('B6A', B6Als8);

var B7Mls8 = ee.Image(post2.get('RADIANCE_MULT_BAND_7'));
var B7Als8 = ee.Image(post2.get('RADIANCE_ADD_BAND_7'));
print ('B7M', B7Mls8);
print ('B7A', B7Als8);

var BMls8 = [B2Mls8, B3Mls8, B4Mls8, B5Mls8, B6Mls8, B7Mls8];
var BAls8 = [B2Als8, B3Als8, B4Als8, B5Als8, B6Als8, B7Als8];
print (BMls8);
print (BAls8);

//create an empty image to save new radiance bands to
var ls8toaRad = ee.Image(0);

//calculate radiance by looping through bands
for (var i=0;i<bandsls8.length; i++){
  
  var band = post.select(bandsls8[i]);
  print (band);
  var mult = ee.Number(BMls8[i]);
  print(mult);
  var add = ee.Number(BAls8[i]);
  print(add);
  var Radls8 = (band.multiply(ee.Image(mult)).add(ee.Image(add)));
  
  ls8toaRad = ls8toaRad.addBands(Radls8.rename(bandsls8[i]));
}  

//check result
print (ls8toaRad);

//display result
var trueColor2 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 100};
Map.addLayer(ls8toaRad, trueColor2,'Post Landsat Radiance True Color');  
  
////////////////////////////////////////////////////////////////////////////////////////
//define function for atmospheric correction to go from TOA Radiance to surface reflectance for Landsat 5

function costzls5(img) {
  // implements the methods from http://info.asprs.org/publications/pers/96journal/september/1996_sep_1025-1036.pdf
  
  var Esun = {'TM':[1958, 1827, 1551, 1036, 214.9, 80.65],};//Esun values taken from: https://landsat.usgs.gov/esun
  
  var sunAngle = pre.metadata('SUN_ELEVATION','sunAngle').subtract(ee.Image(90));
  var esunVal = Esun['TM'];
  
  var bands = ['B1', 'B2','B3','B4','B5','B7'];
  var obsDate = ee.Date(pre.get('system:time_start'));

  var jan01 = ee.Date.fromYMD(obsDate.get('year'),1,1);
  var doy_index = obsDate.difference(jan01,'day').toUint16();
  
  var dsun2 = ee.Image(1).subtract(ee.Image(0.01672).multiply(ee.Image(0.9856).multiply(ee.Image(doy_index).subtract(ee.Image(4)))).cos()).pow(ee.Image(2));
  
  var Presr = ee.Image(0);
  var piImg = ee.Image(Math.PI);
  
  for (var i=0;i<bands.length;i++){
    var esunImg = ee.Image(esunVal[i]);

    var band = img.select(bands[i]);
    //var denom = esunImg.multiply(sunAngle.multiply(piImg.divide(180)).cos());
    //var numer = piImg.multiply(dsun2);
    var radiance = band.rename('rad');
    
    var darkObject = radiance.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:post.geometry(),scale:150,bestEffort:true});
    var doImg = ee.Image(ee.Number(darkObject.get('rad')));

    var temp = ee.Image(0.01).multiply(esunImg).multiply(sunAngle.multiply(piImg.divide(180)).cos()).divide(dsun2.multiply(piImg));
    var Lhaze = doImg.subtract(temp);

    var band_sr = piImg.multiply(band.subtract(Lhaze)).multiply(dsun2).divide(esunImg.multiply(sunAngle.multiply(piImg.divide(180)).cos()));
    
    Presr = Presr.addBands(band_sr.rename(bands[i]));
  }
  
  return Presr.set('system:time_start',obsDate);
}



////////////////////////////////////////////////////////////////////////////////////////
//define function for atmospheric correction to go from TOA Radiance to surface reflectance for Landsat 8

function costzls8(img) {
  // implements the methods from http://info.asprs.org/publications/pers/96journal/september/1996_sep_1025-1036.pdf
  
  var Esun = {'OLI':[2004.57,1820.75,1549.49, 951.76,247.55,85.46],};

  var sunAngle = post.metadata('SUN_ELEVATION','sunAngle').subtract(ee.Image(90));
  var esunVal = Esun['OLI'];
  
  var bands = ['B2','B3','B4','B5','B6','B7'];
  var obsDate = ee.Date(post.get('system:time_start'));

  var jan01 = ee.Date.fromYMD(obsDate.get('year'),1,1);
  var doy_index = obsDate.difference(jan01,'day').toUint16();
  
  var dsun2 = ee.Image(1).subtract(ee.Image(0.01672).multiply(ee.Image(0.9856).multiply(ee.Image(doy_index).subtract(ee.Image(4)))).cos()).pow(ee.Image(2));
  
  var Postsr = ee.Image(0);
  var piImg = ee.Image(Math.PI);
  
  for (var i=0;i<bands.length;i++){
    var esunImg = ee.Image(esunVal[i]);

    var band = img.select(bands[i]);
    //var denom = esunImg.multiply(sunAngle.multiply(piImg.divide(180)).cos());
    //var numer = piImg.multiply(dsun2);
    var radiance = band.rename('rad');
    
    //change radiance to band
    var darkObject = radiance.reduceRegion({reducer:ee.Reducer.min(),geometry:post.geometry(),scale:150,bestEffort:true});
    var doImg = ee.Image(ee.Number(darkObject.get('rad')));

    var temp = ee.Image(0.01).multiply(esunImg).multiply(sunAngle.multiply(piImg.divide(180)).cos()).divide(dsun2.multiply(piImg));
    var Lhaze = doImg.subtract(temp);

    //change radiance to band
    var band_sr = piImg.multiply(band.subtract(Lhaze)).multiply(dsun2).divide(esunImg.multiply(sunAngle.multiply(piImg.divide(180)).cos()));
    
    Postsr = Postsr.addBands(band_sr.rename(bands[i]));
  }
  
  return Postsr.set('system:time_start',obsDate);
}


//*****************APPLY atmospheric correction****************////////////////////
var preCOST = costzls5(ls5toaRad);
var postCOST = costzls8(ls8toaRad);

var trueColor10 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.5};
var trueColor11 = {bands: ['B3', 'B2', 'B1'], min: 0, max: 0.5};
Map.addLayer(preCOST,trueColor11,'pre COST corrected');
Map.addLayer(postCOST,trueColor10,'post COST corrected');

// ************** APPLY BRDF CORRECTION ***************

var preBRDF = brdfCorrectls5(preCOST);
var postBRDF = brdfCorrectls8(postCOST);

Map.addLayer(preBRDF, trueColor11, 'Pre BRDF');
Map.addLayer(postBRDF, trueColor10, 'Post BRDF');
  
// ************** BRDF CORRECTION L5 **************
  
    function brdfCorrectls5(image) {
      var inputBandNames = image.bandNames(); 
      var constants = {
        pi: Math.PI
      };
      print (inputBandNames);
      var coefficientsByBand = {
        'B1': {fiso: 0.0774, fgeo: 0.0079, fvol: 0.0372},
        'B2': {fiso: 0.1306, fgeo: 0.0178, fvol: 0.0580},
        'B3': {fiso: 0.1690, fgeo: 0.0227, fvol: 0.0574},
        'B4': {fiso: 0.3093, fgeo: 0.0330, fvol: 0.1535},
        'B5': {fiso: 0.3430, fgeo: 0.0453, fvol: 0.1154},
        'B7': {fiso: 0.2658, fgeo: 0.0387, fvol: 0.0639}
      };
      
      var corners = findCorners();
      print (corners);
      
     
      
      viewAngles();
      solarPosition();
      sunZenOut();
      set('relativeSunViewAz', 'i.sunAz - i.viewAz');
      rossThick('kvol', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz');
      rossThick('kvol0', 'i.sunZenOut', 0, 0);
      liThin('kgeo', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz');
      liThin('kgeo0', 'i.sunZenOut', 0, 0);
      adjustBands();
      return image.select(inputBandNames);
      print(inputBandNames);

      function viewAngles() {
        var maxDistanceToSceneEdge = 1000000;
        var maxSatelliteZenith = 7.5;
        var upperCenter = pointBetween(corners.upperLeft, corners.upperRight);
        var lowerCenter = pointBetween(corners.lowerLeft, corners.lowerRight);
        var slope = slopeBetween(lowerCenter, upperCenter);
        var slopePerp = ee.Number(-1).divide(slope);
        set('viewAz',
          ee.Image(ee.Number(Math.PI / 2).subtract((slopePerp).atan())));
    
        var leftLine = toLine(corners.upperLeft, corners.lowerLeft);
        var rightLine = toLine(corners.upperRight, corners.lowerRight);
        var leftDistance = ee.FeatureCollection(leftLine).distance(maxDistanceToSceneEdge);
        var rightDistance = ee.FeatureCollection(rightLine).distance(maxDistanceToSceneEdge);
        var viewZenith = rightDistance.multiply(maxSatelliteZenith * 2)
          .divide(rightDistance.add(leftDistance))
          .subtract(maxSatelliteZenith);
        set('viewZen',
          viewZenith.multiply(Math.PI).divide(180));
      }
        print (viewZenith);
    
      function solarPosition() {
        // Ported from http://pythonfmask.org/en/latest/_modules/fmask/landsatangles.html
        var date = ee.Date(ee.Number(pre.get('system:time_start')));
        var secondsInHour = 3600;
        set('longDeg',
          ee.Image.pixelLonLat().select('longitude'));
        set('latRad',
          ee.Image.pixelLonLat().select('latitude')
            .multiply(Math.PI).divide(180));
        set('hourGMT',
          ee.Number(date.getRelative('second', 'day')).divide(secondsInHour));
        set('jdp', // Julian Date Proportion
          date.getFraction('year'));
        set('jdpr', // Julian Date Proportion in Radians
          'i.jdp * 2 * {pi}');
        set('meanSolarTime',
          'i.hourGMT + i.longDeg / 15');
        set('localSolarDiff',
          '(0.000075 + 0.001868 * cos(i.jdpr) - 0.032077 * sin(i.jdpr)' +
          '- 0.014615 * cos(2 * i.jdpr) - 0.040849 * sin(2 * i.jdpr))' +
          '* 12 * 60 / {pi}');
        set('trueSolarTime',
          'i.meanSolarTime + i.localSolarDiff / 60 - 12');
        set('angleHour',
          'i.trueSolarTime * 15 * {pi} / 180');
        set('delta',
          '0.006918 - 0.399912 * cos(i.jdpr) + 0.070257 * sin(i.jdpr) - 0.006758 * cos(2 * i.jdpr)' +
          '+ 0.000907 * sin(2 * i.jdpr) - 0.002697 * cos(3 * i.jdpr) + 0.001480 * sin(3 * i.jdpr)');
        set('cosSunZen',
          'sin(i.latRad) * sin(i.delta) ' +
          '+ cos(i.latRad) * cos(i.delta) * cos(i.angleHour)');
        set('sunZen',
          'acos(i.cosSunZen)');
        set('sinSunAzSW',
          toImage('cos(i.delta) * sin(i.angleHour) / sin(i.sunZen)')
            .clamp(-1, 1));
        set('cosSunAzSW',
          '(-cos(i.latRad) * sin(i.delta)' +
          '+ sin(i.latRad) * cos(i.delta) * cos(i.angleHour)) / sin(i.sunZen)');
        set('sunAzSW',
          'asin(i.sinSunAzSW)');
        setIf('sunAzSW',
          'i.cosSunAzSW <= 0',
          '{pi} - i.sunAzSW',
          'sunAzSW');
        setIf('sunAzSW',
          'i.cosSunAzSW > 0 and i.sinSunAzSW <= 0',
          '2 * {pi} + i.sunAzSW',
          'sunAzSW');
        set('sunAz',
          'i.sunAzSW + {pi}');
        setIf('sunAz',
          'i.sunAz > 2 * {pi}',
          'i.sunAz - 2 * {pi}',
          'sunAz');
      }
    
      function sunZenOut() {
        // https://nex.nasa.gov/nex/static/media/publication/HLS.v1.0.UserGuide.pdf
        set('centerLat',
          ee.Number(
            ee.Geometry(pre.get('system:footprint'))
              .bounds().centroid(30).coordinates().get(0))
            .multiply(Math.PI).divide(180));
        set('sunZenOut',
          '(31.0076' +
          '- 0.1272 * i.centerLat' +
          '+ 0.01187 * pow(i.centerLat, 2)' +
          '+ 2.40E-05 * pow(i.centerLat, 3)' +
          '- 9.48E-07 * pow(i.centerLat, 4)' +
          '- 1.95E-09 * pow(i.centerLat, 5)' +
          '+ 6.15E-11 * pow(i.centerLat, 6)) * {pi}/180');
      }
    
      function rossThick(bandName, sunZen, viewZen, relativeSunViewAz) {
        var args = {sunZen: sunZen, viewZen: viewZen, relativeSunViewAz: relativeSunViewAz};
        cosPhaseAngle('cosPhaseAngle', sunZen, viewZen, relativeSunViewAz);
        set('phaseAngle',
          'acos(i.cosPhaseAngle)');
        set(bandName,
          '(({pi}/2 - i.phaseAngle) * i.cosPhaseAngle + sin(i.phaseAngle)) ' +
          '/ (cos({sunZen}) + cos({viewZen})) - {pi}/4', args);
      }
    
      function liThin(bandName, sunZen, viewZen, relativeSunViewAz) {
        // From https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf
        var args = {
          sunZen: sunZen,
          viewZen: viewZen,
          relativeSunViewAz: relativeSunViewAz,
          'h/b': 2,
        };
    
        anglePrime('sunZenPrime', sunZen);
        anglePrime('viewZenPrime', viewZen);
        cosPhaseAngle('cosPhaseAnglePrime', 'i.sunZenPrime', 'i.viewZenPrime', relativeSunViewAz);
        set('distance',
          'sqrt(pow(tan(i.sunZenPrime), 2) + pow(tan(i.viewZenPrime), 2)' +
          '- 2 * tan(i.sunZenPrime) * tan(i.viewZenPrime) * cos({relativeSunViewAz}))', args);
        set('temp',
          '1/cos(i.sunZenPrime) + 1/cos(i.viewZenPrime)');
        set('cosT',
          toImage('{h/b} * sqrt(pow(i.distance, 2) + pow(tan(i.sunZenPrime) * tan(i.viewZenPrime) * sin({relativeSunViewAz}), 2))' +
            '/ i.temp', args)
            .clamp(-1, 1));
        set('t', 'acos(i.cosT)');
        set('overlap',
          '(1/{pi}) * (i.t - sin(i.t) * i.cosT) * (i.temp)');
        setIf('overlap', 'i.overlap > 0', 0);
        set(bandName,
          'i.overlap - i.temp' +
          '+ (1/2) * (1 + i.cosPhaseAnglePrime) * (1/cos(i.sunZenPrime)) * (1/cos(i.viewZenPrime))');
      }
    
      function anglePrime(name, angle) {
        var args = {'b/r': 1, angle: angle};
        set('tanAnglePrime',
          '{b/r} * tan({angle})', args);
        setIf('tanAnglePrime', 'i.tanAnglePrime < 0', 0);
        set(name,
          'atan(i.tanAnglePrime)');
      }
    
      function cosPhaseAngle(name, sunZen, viewZen, relativeSunViewAz) {
        var args = {
          sunZen: sunZen,
          viewZen: viewZen,
          relativeSunViewAz: relativeSunViewAz
        };
        set(name,
          toImage('cos({sunZen}) * cos({viewZen})' +
            '+ sin({sunZen}) * sin({viewZen}) * cos({relativeSunViewAz})', args)
            .clamp(-1, 1));
      }
    
      function adjustBands() {
        for (var bandName in coefficientsByBand)
          applyCFactor(bandName, coefficientsByBand[bandName]);
      }
    
      function applyCFactor(bandName, coefficients) {
        brdf('brdf', 'kvol', 'kgeo', coefficients);
        brdf('brdf0', 'kvol0', 'kgeo0', coefficients);
        set('cFactor',
          'i.brdf0 / i.brdf', coefficients);
        set(bandName,
          '{bandName} * i.cFactor', {bandName: 'i.' + bandName});
      }
    
      function brdf(bandName, kvolBand, kgeoBand, coefficients) {
        var args = merge(coefficients, {
          // kvol: 'i.' + kvolBand,
          kvol: '3 * i.' + kvolBand,     // check this multiplication factor.  Is there an 'optimal' value?  Without a factor here, there is not enough correction.
          kgeo: 'i.' + kgeoBand
        });
        return set(bandName,
          '{fiso} + {fvol} * {kvol} + {fgeo} * {kvol}', args);
      }
    
      function findCorners() {
        var footprint = ee.Geometry(pre.get('system:footprint'));
        var bounds = ee.List(footprint.bounds().coordinates().get(0));
        var coords = footprint.coordinates();
    
        var xs = coords.map(function (item) {
          return x(item);
        });
        var ys = coords.map(function (item) {
          return y(item);
        });
    
        function findCorner(targetValue, values) {
          var diff = values.map(function (value) {
            return ee.Number(value).subtract(targetValue).abs();
          });
          var minValue = diff.reduce(ee.Reducer.min());
          var idx = diff.indexOf(minValue);
          return coords.get(idx);
        }
    
        var lowerLeft = findCorner(x(bounds.get(0)), xs);
        var lowerRight = findCorner(y(bounds.get(1)), ys);
        var upperRight = findCorner(x(bounds.get(2)), xs);
        var upperLeft = findCorner(y(bounds.get(3)), ys);
        return {
          upperLeft: upperLeft,
          upperRight: upperRight,
          lowerRight: lowerRight,
          lowerLeft: lowerLeft
        };
      }
  
      function x(point) {
        return ee.Number(ee.List(point).get(0));
      }
    
      function y(point) {
        return ee.Number(ee.List(point).get(1));
      }
    
      function pointBetween(pointA, pointB) {
        return ee.Geometry.LineString([pointA, pointB]).centroid().coordinates();
      }
    
      function slopeBetween(pointA, pointB) {
        return ((y(pointA)).subtract(y(pointB))).divide((x(pointA)).subtract(x(pointB)));
      }
    
      function toLine(pointA, pointB) {
        return ee.Geometry.LineString([pointA, pointB]);
      }

// ************** COMMON HELPERS **************

      function set(name, toAdd, args) {
        toAdd = toImage(toAdd, args)
        image = image.addBands(toAdd.rename(name), null, true)
      }
    
      function setIf(name, condition, trueValue, falseValue) {
        condition = toImage(condition)
        var trueMasked = toImage(trueValue).mask(toImage(condition))
        var falseMasked = toImage(falseValue).mask(invertMask(condition))
        var value = trueMasked.unmask(falseMasked)
        set(name, value)
        
    
        function invertMask(mask) {
          return mask.multiply(-1).add(1)
        }
      }
    
      function toImage(band, args) {
        if ((typeof band) === 'string') {
          if (band.indexOf('.') > -1 || band.indexOf(' ') > -1 || band.indexOf('{') > -1) {
            band = image.expression(format(band, args), {i: image})
          } else
            band = image.select(band)
        }
        return ee.Image(band)
      }
    
      function format(s, args) {
        if (!args) args = {}
        var allArgs = merge(constants, args)
        var result = s.replace(/{([^{}]*)}/g,
          function (a, b) {
            var replacement = allArgs[b]
            if (replacement == null) {
              print('Undeclared argument: ' + b, 's: ' + s, args)
              return null
            }
            return allArgs[b]
          }
        )
        if (result.indexOf('{') > -1)
          return format(result, args)
        return result
      }
      
      function merge(o1, o2) {
        function addAll(target, toAdd) {
          for (var key in toAdd) target[key] = toAdd[key]
        }
    
        var result = {}
        addAll(result, o1)
        addAll(result, o2)
        return result
      }
    
      function show(band, min, max) {
        Map.addLayer(toImage(band), {min: min ? min : -1, max: max ? max : 1}, band)
      }
    }
    

// ************** BRDF CORRECTION L8 **************
  
    function brdfCorrectls8(image) {
      var inputBandNames = image.bandNames() 
      var constants = {
        pi: Math.PI
      }
      print (inputBandNames)
      var coefficientsByBand = {
        'B2': {fiso: 0.0774, fgeo: 0.0079, fvol: 0.0372},
        'B3': {fiso: 0.1306, fgeo: 0.0178, fvol: 0.0580},
        'B4': {fiso: 0.1690, fgeo: 0.0227, fvol: 0.0574},
        'B5': {fiso: 0.3093, fgeo: 0.0330, fvol: 0.1535},
        'B6': {fiso: 0.3430, fgeo: 0.0453, fvol: 0.1154},
        'B7': {fiso: 0.2658, fgeo: 0.0387, fvol: 0.0639}
      }
      
      var corners = findCorners()
      print (corners)
      
     
      
      viewAngles()
      solarPosition()
      sunZenOut()
      set('relativeSunViewAz', 'i.sunAz - i.viewAz')
      rossThick('kvol', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz')
      rossThick('kvol0', 'i.sunZenOut', 0, 0)
      liThin('kgeo', 'i.sunZen', 'i.viewZen', 'i.relativeSunViewAz')
      liThin('kgeo0', 'i.sunZenOut', 0, 0)
      adjustBands()
      return image.select(inputBandNames)
      print(inputBandNames)

      function viewAngles() {
        var maxDistanceToSceneEdge = 1000000
        var maxSatelliteZenith = 7.5
        var upperCenter = pointBetween(corners.upperLeft, corners.upperRight)
        var lowerCenter = pointBetween(corners.lowerLeft, corners.lowerRight)
        var slope = slopeBetween(lowerCenter, upperCenter)
        var slopePerp = ee.Number(-1).divide(slope)
        set('viewAz',
          ee.Image(ee.Number(Math.PI / 2).subtract((slopePerp).atan())))
    
        var leftLine = toLine(corners.upperLeft, corners.lowerLeft)
        var rightLine = toLine(corners.upperRight, corners.lowerRight)
        var leftDistance = ee.FeatureCollection(leftLine).distance(maxDistanceToSceneEdge)
        var rightDistance = ee.FeatureCollection(rightLine).distance(maxDistanceToSceneEdge)
        var viewZenith = rightDistance.multiply(maxSatelliteZenith * 2)
          .divide(rightDistance.add(leftDistance))
          .subtract(maxSatelliteZenith)
        set('viewZen',
          viewZenith.multiply(Math.PI).divide(180))
      }
        print (viewZenith)
    
      function solarPosition() {
        // Ported from http://pythonfmask.org/en/latest/_modules/fmask/landsatangles.html
        var date = ee.Date(ee.Number(post.get('system:time_start')))
        var secondsInHour = 3600
        set('longDeg',
          ee.Image.pixelLonLat().select('longitude'))
        set('latRad',
          ee.Image.pixelLonLat().select('latitude')
            .multiply(Math.PI).divide(180))
        set('hourGMT',
          ee.Number(date.getRelative('second', 'day')).divide(secondsInHour))
        set('jdp', // Julian Date Proportion
          date.getFraction('year'))
        set('jdpr', // Julian Date Proportion in Radians
          'i.jdp * 2 * {pi}')
        set('meanSolarTime',
          'i.hourGMT + i.longDeg / 15')
        set('localSolarDiff',
          '(0.000075 + 0.001868 * cos(i.jdpr) - 0.032077 * sin(i.jdpr)' +
          '- 0.014615 * cos(2 * i.jdpr) - 0.040849 * sin(2 * i.jdpr))' +
          '* 12 * 60 / {pi}')
        set('trueSolarTime',
          'i.meanSolarTime + i.localSolarDiff / 60 - 12')
        set('angleHour',
          'i.trueSolarTime * 15 * {pi} / 180')
        set('delta',
          '0.006918 - 0.399912 * cos(i.jdpr) + 0.070257 * sin(i.jdpr) - 0.006758 * cos(2 * i.jdpr)' +
          '+ 0.000907 * sin(2 * i.jdpr) - 0.002697 * cos(3 * i.jdpr) + 0.001480 * sin(3 * i.jdpr)')
        set('cosSunZen',
          'sin(i.latRad) * sin(i.delta) ' +
          '+ cos(i.latRad) * cos(i.delta) * cos(i.angleHour)')
        set('sunZen',
          'acos(i.cosSunZen)')
        set('sinSunAzSW',
          toImage('cos(i.delta) * sin(i.angleHour) / sin(i.sunZen)')
            .clamp(-1, 1))
        set('cosSunAzSW',
          '(-cos(i.latRad) * sin(i.delta)' +
          '+ sin(i.latRad) * cos(i.delta) * cos(i.angleHour)) / sin(i.sunZen)')
        set('sunAzSW',
          'asin(i.sinSunAzSW)')
        setIf('sunAzSW',
          'i.cosSunAzSW <= 0',
          '{pi} - i.sunAzSW',
          'sunAzSW')
        setIf('sunAzSW',
          'i.cosSunAzSW > 0 and i.sinSunAzSW <= 0',
          '2 * {pi} + i.sunAzSW',
          'sunAzSW')
        set('sunAz',
          'i.sunAzSW + {pi}')
        setIf('sunAz',
          'i.sunAz > 2 * {pi}',
          'i.sunAz - 2 * {pi}',
          'sunAz')
      }
    
      function sunZenOut() {
        // https://nex.nasa.gov/nex/static/media/publication/HLS.v1.0.UserGuide.pdf
        set('centerLat',
          ee.Number(
            ee.Geometry(post.get('system:footprint'))
              .bounds().centroid(30).coordinates().get(0))
            .multiply(Math.PI).divide(180))
        set('sunZenOut',
          '(31.0076' +
          '- 0.1272 * i.centerLat' +
          '+ 0.01187 * pow(i.centerLat, 2)' +
          '+ 2.40E-05 * pow(i.centerLat, 3)' +
          '- 9.48E-07 * pow(i.centerLat, 4)' +
          '- 1.95E-09 * pow(i.centerLat, 5)' +
          '+ 6.15E-11 * pow(i.centerLat, 6)) * {pi}/180')
      }
    
      function rossThick(bandName, sunZen, viewZen, relativeSunViewAz) {
        var args = {sunZen: sunZen, viewZen: viewZen, relativeSunViewAz: relativeSunViewAz}
        cosPhaseAngle('cosPhaseAngle', sunZen, viewZen, relativeSunViewAz)
        set('phaseAngle',
          'acos(i.cosPhaseAngle)')
        set(bandName,
          '(({pi}/2 - i.phaseAngle) * i.cosPhaseAngle + sin(i.phaseAngle)) ' +
          '/ (cos({sunZen}) + cos({viewZen})) - {pi}/4', args)
      }
    
      function liThin(bandName, sunZen, viewZen, relativeSunViewAz) {
        // From https://modis.gsfc.nasa.gov/data/atbd/atbd_mod09.pdf
        var args = {
          sunZen: sunZen,
          viewZen: viewZen,
          relativeSunViewAz: relativeSunViewAz,
          'h/b': 2,
        }
    
        anglePrime('sunZenPrime', sunZen)
        anglePrime('viewZenPrime', viewZen)
        cosPhaseAngle('cosPhaseAnglePrime', 'i.sunZenPrime', 'i.viewZenPrime', relativeSunViewAz)
        set('distance',
          'sqrt(pow(tan(i.sunZenPrime), 2) + pow(tan(i.viewZenPrime), 2)' +
          '- 2 * tan(i.sunZenPrime) * tan(i.viewZenPrime) * cos({relativeSunViewAz}))', args)
        set('temp',
          '1/cos(i.sunZenPrime) + 1/cos(i.viewZenPrime)')
        set('cosT',
          toImage('{h/b} * sqrt(pow(i.distance, 2) + pow(tan(i.sunZenPrime) * tan(i.viewZenPrime) * sin({relativeSunViewAz}), 2))' +
            '/ i.temp', args)
            .clamp(-1, 1))
        set('t', 'acos(i.cosT)')
        set('overlap',
          '(1/{pi}) * (i.t - sin(i.t) * i.cosT) * (i.temp)')
        setIf('overlap', 'i.overlap > 0', 0)
        set(bandName,
          'i.overlap - i.temp' +
          '+ (1/2) * (1 + i.cosPhaseAnglePrime) * (1/cos(i.sunZenPrime)) * (1/cos(i.viewZenPrime))')
      }
    
      function anglePrime(name, angle) {
        var args = {'b/r': 1, angle: angle}
        set('tanAnglePrime',
          '{b/r} * tan({angle})', args)
        setIf('tanAnglePrime', 'i.tanAnglePrime < 0', 0)
        set(name,
          'atan(i.tanAnglePrime)')
      }
    
      function cosPhaseAngle(name, sunZen, viewZen, relativeSunViewAz) {
        var args = {
          sunZen: sunZen,
          viewZen: viewZen,
          relativeSunViewAz: relativeSunViewAz
        }
        set(name,
          toImage('cos({sunZen}) * cos({viewZen})' +
            '+ sin({sunZen}) * sin({viewZen}) * cos({relativeSunViewAz})', args)
            .clamp(-1, 1))
      }
    
      function adjustBands() {
        for (var bandName in coefficientsByBand)
          applyCFactor(bandName, coefficientsByBand[bandName])
      }
    
      function applyCFactor(bandName, coefficients) {
        brdf('brdf', 'kvol', 'kgeo', coefficients)
        brdf('brdf0', 'kvol0', 'kgeo0', coefficients)
        set('cFactor',
          'i.brdf0 / i.brdf', coefficients)
        set(bandName,
          '{bandName} * i.cFactor', {bandName: 'i.' + bandName})
      }
    
      function brdf(bandName, kvolBand, kgeoBand, coefficients) {
        var args = merge(coefficients, {
          // kvol: 'i.' + kvolBand,
          kvol: '3 * i.' + kvolBand,     // check this multiplication factor.  Is there an 'optimal' value?  Without a factor here, there is not enough correction.
          kgeo: 'i.' + kgeoBand
        })
        return set(bandName,
          '{fiso} + {fvol} * {kvol} + {fgeo} * {kvol}', args)
      }
    
      function findCorners() {
        var footprint = ee.Geometry(post.get('system:footprint'))
        var bounds = ee.List(footprint.bounds().coordinates().get(0))
        var coords = footprint.coordinates()
    
        var xs = coords.map(function (item) {
          return x(item)
        })
        var ys = coords.map(function (item) {
          return y(item);
        });
    
        function findCorner(targetValue, values) {
          var diff = values.map(function (value) {
            return ee.Number(value).subtract(targetValue).abs();
          });
          var minValue = diff.reduce(ee.Reducer.min());
          var idx = diff.indexOf(minValue);
          return coords.get(idx);
        }
    
        var lowerLeft = findCorner(x(bounds.get(0)), xs);
        var lowerRight = findCorner(y(bounds.get(1)), ys);
        var upperRight = findCorner(x(bounds.get(2)), xs);
        var upperLeft = findCorner(y(bounds.get(3)), ys);
        return {
          upperLeft: upperLeft,
          upperRight: upperRight,
          lowerRight: lowerRight,
          lowerLeft: lowerLeft
        };
      }
  
      function x(point) {
        return ee.Number(ee.List(point).get(0));
      }
    
      function y(point) {
        return ee.Number(ee.List(point).get(1));
      }
    
      function pointBetween(pointA, pointB) {
        return ee.Geometry.LineString([pointA, pointB]).centroid().coordinates();
      }
    
      function slopeBetween(pointA, pointB) {
        return ((y(pointA)).subtract(y(pointB))).divide((x(pointA)).subtract(x(pointB)));
      }
    
      function toLine(pointA, pointB) {
        return ee.Geometry.LineString([pointA, pointB]);
      }

// ************** COMMON HELPERS **************

      function set(name, toAdd, args) {
        toAdd = toImage(toAdd, args)
        image = image.addBands(toAdd.rename(name), null, true)
      }
    
      function setIf(name, condition, trueValue, falseValue) {
        condition = toImage(condition)
        var trueMasked = toImage(trueValue).mask(toImage(condition))
        var falseMasked = toImage(falseValue).mask(invertMask(condition))
        var value = trueMasked.unmask(falseMasked)
        set(name, value)
        
    
        function invertMask(mask) {
          return mask.multiply(-1).add(1)
        }
      }
    
      function toImage(band, args) {
        if ((typeof band) === 'string') {
          if (band.indexOf('.') > -1 || band.indexOf(' ') > -1 || band.indexOf('{') > -1) {
            band = image.expression(format(band, args), {i: image})
          } else
            band = image.select(band)
        }
        return ee.Image(band)
      }
    
      function format(s, args) {
        if (!args) args = {}
        var allArgs = merge(constants, args)
        var result = s.replace(/{([^{}]*)}/g,
          function (a, b) {
            var replacement = allArgs[b]
            if (replacement == null) {
              print('Undeclared argument: ' + b, 's: ' + s, args)
              return null
            }
            return allArgs[b]
          }
        )
        if (result.indexOf('{') > -1)
          return format(result, args)
        return result
      }
      
      function merge(o1, o2) {
        function addAll(target, toAdd) {
          for (var key in toAdd) target[key] = toAdd[key]
        }
    
        var result = {}
        addAll(result, o1)
        addAll(result, o2)
        return result
      }
    
      function show(band, min, max) {
        Map.addLayer(toImage(band), {min: min ? min : -1, max: max ? max : 1}, band)
      }
    }
    

/////////////////////////////////////////////////////////////////////    
//var dif = postBRDF.subtract(postCOST)
//Map.addLayer (dif)
/////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////////////////////////
//Get values of spectral endmemebers - these are taken from pre-image (Landsat 5) after atmospheric correction and BRDF

var soilDict = postBRDF.reduceRegion(ee.Reducer.mean(),soilL5ACl8, 30);
var soilList = soilDict.values(['B2', 'B3', 'B4','B5', 'B6', 'B7']);
print (soilDict);
print (soilList);

var vegDict = postBRDF.reduceRegion(ee.Reducer.mean(),vegL5AC, 30);
var vegList = vegDict.values(['B2', 'B3', 'B4','B5', 'B6', 'B7']);
print (vegDict);
print (vegList);

var waterDict = postBRDF.reduceRegion(ee.Reducer.mean(),waterL5AC, 30);
var waterList = waterDict.values(['B2', 'B3', 'B4','B5', 'B6', 'B7']);
print (waterDict);
print (waterList);

///////////////////////////////////////////////////////////////////////
//select the post bands to be used in spectral unmixing
var postACb = postBRDF.select('B2', 'B3','B4', 'B5', 'B6','B7');
// Unmix post image.
var fractionsPost = postACb.unmix([soilList, vegList, waterList], true,true);
Map.addLayer(fractionsPost, {min:0, max: 0.6}, 'Post Unmixed Image');

//select pre bands for spectral unmixing
var preACb = preBRDF.select('B1', 'B2', 'B3','B4', 'B5', 'B7');
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
var mndwiPre = preBRDF.normalizedDifference(['B2', 'B5']);
var mndwiPost = postBRDF.normalizedDifference(['B3','B6']);

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
//display "gage" values
var geometry2 = [gage1, gage3, gage4, gage5, gage6, gage7, gage10, gage11, gage12, gage13, gage14, gage16, gage18, gage19, gage20, gage21, gage22, gage23, gage25];
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
