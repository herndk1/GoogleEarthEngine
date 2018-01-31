//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//Preliminary: evaluating the impact of various methods of atmospheric correction on the spectral signatures of various land covers
//Kelsey E. Herndon, Graduate Research Assistant, SERVIR Science Coordination Office, The University of Alabama in Huntsville
//Date created: 9/15/2017
//Last modified: 11/1/2017 
//Purpose: This script evaluates the impact of different atmospheric correction methods on detecting surface water in the Ferlo Region of Senegal. Landsat TOA 
//reflectance, Landsat USGS surface reflectance (6S), and surface reflectance calculated using the COST method are assessed. 
//Notes: Kel Markert (NASA SERVIR) provided the COST atmospheric correction function that is used in this script.
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

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
          [-16.53167724609375, 16.45715879614139]]]),
    Water = /* color: #bf04c2 */ee.Geometry.Polygon(
        [[[-15.822372436523438, 16.325411207783855],
          [-15.82305908203125, 16.315526551275934],
          [-15.810012817382812, 16.314867556423422],
          [-15.810012817382812, 16.32672912425378]]]),
    Desert = /* color: #ff0000 */ee.Geometry.Polygon(
        [[[-15.433731079101562, 16.378779715441343],
          [-15.433731079101562, 16.368897759211915],
          [-15.411758422851562, 16.366921307880315],
          [-15.413818359375, 16.38009727176448]]]),
    Veg = /* color: #00ff00 */ee.Geometry.Polygon(
        [[[-15.68796157836914, 16.40027124032415],
          [-15.689420700073242, 16.397389372748126],
          [-15.684871673583984, 16.39516618861569],
          [-15.683069229125977, 16.397636391640415]]]),
    Coast = /* color: #0000ff */ee.Geometry.Polygon(
        [[[-16.496872901916504, 16.347402776150563],
          [-16.499361991882324, 16.34662034707395],
          [-16.500906944274902, 16.3408961653517],
          [-16.49867534637451, 16.34176098053851]]]);

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

////////////////////////////////////////////////////////////////////////////////////
//define function for atmospheric correction 
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
////////////////////////////////////////////////////////////////////////////////////

//apply COST function to Landsat TOA scene
var lsCOST = costz(lstoa);

//////////////////////////////////////////////////////////////////////////////////// 
//display COST results 
var trueColor1 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3};
var trueColor2 = {bands: ['b4', 'b3', 'b2'], min: 0, max: 3000};
//display COST atmospherically corrected scene (true color)
Map.addLayer(lsCOST, trueColor1, 'Landsat with COST atmospheric correction');

//////////////////////////////////////////////////////////////////////////
//calculate difference between COST and SR for each band
//blue band 
var blueComp = lsCOST.select('B2').subtract(lssr.select('b2').multiply(0.0001));
Map.addLayer(blueComp,{min:-0.5, max: 0.5}, 'Blue Band Comparison');
//green band
var greenComp = lsCOST.select('B3').subtract(lssr.select('b3').multiply(0.0001));
Map.addLayer(greenComp,{min:-0.5, max: 0.5}, 'Green Band Comparison');
//red band
var redComp = lsCOST.select('B4').subtract(lssr.select('b4').multiply(0.0001));
Map.addLayer(redComp,{min:-0.5, max: 0.5}, 'Red Band Comparison');
//nir band
var nirComp = lsCOST.select('B5').subtract(lssr.select('b5').multiply(0.0001));
Map.addLayer(nirComp,{min:-0.5, max: 0.5}, 'NIR Band Comparison');
//swir1 band
var swir1Comp = lsCOST.select('B6').subtract(lssr.select('b6').multiply(0.0001));
Map.addLayer(swir1Comp,{min:-0.5, max: 0.5}, 'SWIR1 Band Comparison');
//swir2 band
var swir2Comp = lsCOST.select('B7').subtract(lssr.select('b7').multiply(0.0001));
Map.addLayer(swir2Comp,{min:-0.5, max: 0.5}, 'SWIR2 Band Comparison');


/////////////////////////////////////////////////////////////////////////
//display spectral responses of different land covers for each atmospheric correction method
//select bands to collect response from
var lsCOSTselect = ee.Image(lsCOST).select(['B[2-7]']);
var lsTOAselect = ee.Image(lstoa).select(['B[2-7]']);
var lsSRselect = ee.Image(lssr).select(['b[2-7]']).multiply(0.0001);
//define central wavelength of each band
var wavelengths = [0.48,0.56,0.65,0.86,1.61,2.2];

//create one image with all bands from above
var image = lsCOSTselect.addBands(lsTOAselect).addBands(lsSRselect);

// Extract each band mean for each AOI to an array
var costArrayWater = lsCOSTselect.reduceRegion(ee.Reducer.toList(), Water, 30)
                 .toArray(lsCOSTselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var toaArrayWater  = lsTOAselect.reduceRegion(ee.Reducer.toList(), Water, 30)
                 .toArray(lsTOAselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var srArrayWater   = lsSRselect.reduceRegion(ee.Reducer.toList(), Water, 30)
                 .toArray(lsSRselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);

// Create arrays for chart
var arrayWater = ee.Array([costArrayWater.toList(),toaArrayWater.toList(),srArrayWater.toList()]);
var allSpectraWater = arrayWater.slice(1, 0);  // For the Y axis.

// Generate and style the chart.
var chartWater = ui.Chart.array.values(allSpectraWater, 1, wavelengths)
    .setChartType('LineChart')
    .setSeriesNames(['COST', 'TOA', 'ESPA'])
    .setOptions({
      title: 'Spectral Response of Water from Different Atmospheric Corrections',
      vAxes: {
        0: {
          title: 'Refelctance [%]'
        }
      },
      hAxis: {
        title: 'Wavelength [mirometers]'
      },
      interpolateNulls: true,
      pointSize: 3,
      lineWidth: 1,
      series:{
        0: {color: 'red'},
        1: {color: 'gray'},
        2: {color: 'black'}
      }
    });

print(chartWater);

/////////////////////////////////////////////////////////////////////////////////
//Create spectral signature comparison for Vegetation
// Extract band mean for AOI to an array
var costArrayVeg = lsCOSTselect.reduceRegion(ee.Reducer.toList(), Veg, 30)
                 .toArray(lsCOSTselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var toaArrayVeg  = lsTOAselect.reduceRegion(ee.Reducer.toList(), Veg, 30)
                 .toArray(lsTOAselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var srArrayVeg   = lsSRselect.reduceRegion(ee.Reducer.toList(), Veg, 30)
                 .toArray(lsSRselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);

// Create arrays for charting.
var arrayVeg = ee.Array([costArrayVeg.toList(),toaArrayVeg.toList(),srArrayVeg.toList()]);
var allSpectraVeg = arrayVeg.slice(1, 0);  // For the Y axis.

// Generate and style the chart.
var chartVeg = ui.Chart.array.values(allSpectraVeg, 1, wavelengths)
    .setChartType('LineChart')
    .setSeriesNames(['COST', 'TOA', 'ESPA'])
    .setOptions({
      title: 'Spectral Response of Vegetation from Different Atmospheric Corrections',
      vAxes: {
        0: {
          title: 'Refelctance [%]'
        }
      },
      hAxis: {
        title: 'Wavelength [mirometers]'
      },
      interpolateNulls: true,
      pointSize: 3,
      lineWidth: 1,
      series:{
        0: {color: 'red'},
        1: {color: 'gray'},
        2: {color: 'black'}
      }
    });

print(chartVeg);
//////////////////////////////////////////////////////////////////////////////////
//create spectral signature for Desert
// Extract band mean for AOI to an array
var costArrayDesert = lsCOSTselect.reduceRegion(ee.Reducer.toList(), Desert, 30)
                 .toArray(lsCOSTselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var toaArrayDesert  = lsTOAselect.reduceRegion(ee.Reducer.toList(), Desert, 30)
                 .toArray(lsTOAselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var srArrayDesert   = lsSRselect.reduceRegion(ee.Reducer.toList(), Desert, 30)
                 .toArray(lsSRselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);

// Create arrays for charting.
var arrayDesert = ee.Array([costArrayDesert.toList(),toaArrayDesert.toList(),srArrayDesert.toList()]);
var allSpectraDesert = arrayDesert.slice(1, 0);  // For the Y axis.

// Generate and style the chart.
var chartDesert = ui.Chart.array.values(allSpectraDesert, 1, wavelengths)
    .setChartType('LineChart')
    .setSeriesNames(['COST', 'TOA', 'ESPA'])
    .setOptions({
      title: 'Spectral Response of Desert from Different Atmospheric Corrections',
      vAxes: {
        0: {
          title: 'Refelctance [%]'
        }
      },
      hAxis: {
        title: 'Wavelength [mirometers]'
      },
      interpolateNulls: true,
      pointSize: 3,
      lineWidth: 1,
      series:{
        0: {color: 'red'},
        1: {color: 'gray'},
        2: {color: 'black'}
      }
    });

print(chartDesert);
/////////////////////////////////////////////////////////////////////////////////
//create spectral signature for Coast
// Extract band mean for AOI to an array
var costArrayCoast = lsCOSTselect.reduceRegion(ee.Reducer.toList(), Coast, 30)
                 .toArray(lsCOSTselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var toaArrayCoast  = lsTOAselect.reduceRegion(ee.Reducer.toList(), Coast, 30)
                 .toArray(lsTOAselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);
var srArrayCoast   = lsSRselect.reduceRegion(ee.Reducer.toList(), Coast, 30)
                 .toArray(lsSRselect.bandNames(),1).reduce(ee.Reducer.mean(),[1]);

// Create arrays for charting.
var arrayCoast = ee.Array([costArrayCoast.toList(),toaArrayCoast.toList(),srArrayCoast.toList()]);
var allSpectraCoast = arrayCoast.slice(1, 0);  // For the Y axis.

// Generate and style the chart.
var chartCoast = ui.Chart.array.values(allSpectraCoast, 1, wavelengths)
    .setChartType('LineChart')
    .setSeriesNames(['COST', 'TOA', 'ESPA'])
    .setOptions({
      title: 'Spectral Response of Coast from Different Atmospheric Corrections',
      vAxes: {
        0: {
          title: 'Refelctance [%]'
        }
      },
      hAxis: {
        title: 'Wavelength [mirometers]'
      },
      interpolateNulls: true,
      pointSize: 3,
      lineWidth: 1,
      series:{
        0: {color: 'red'},
        1: {color: 'gray'},
        2: {color: 'black'}
      }
    });

print(chartCoast);
