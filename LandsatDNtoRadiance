////////////////////////////////////////////////////////////////
//Kelsey E. Herndon, Graduate Research Assistant, SERVIR Science Coordination Office, The University of Alabama in Huntsville
//kelsey.e.herndon@nasa.gov
//created: 11/8/2017
//last modified: 11/13/2017 
//purpose: To loop through bands of Landsat 8 DN and convert to TOA Radiance 
//using gain/bias in metadata instead of using GEE function
//////////////////////////////////////////////////////////////

//define imports
var l8Raw = ee.ImageCollection("LANDSAT/LC8_L1T");

//define point to filter TOA landsat scenes
var pt = /* color: #d63000 */ee.Geometry.Point([-15.785980224609375, 16.337272136282255]);

//filter landsat scenes
var l8RawSelect = ee.Image(l8Raw
  .filterBounds(pt)
  .filterDate('2016-12-11', '2016-12-13')
  .first());

//display resulting Raw Landsat scene
var trueColor1 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 30000};
Map.addLayer(l8RawSelect, trueColor1,'Landsat Raw true color');

//define bands to be converted to radiance
var bands = ['B1','B2','B3','B4','B5','B6','B7'];


//convert DN to Radiance
//define variables 
var B1M = ee.Image(l8RawSelect.get('RADIANCE_MULT_BAND_1'));
var B1A = ee.Image(l8RawSelect.get('RADIANCE_ADD_BAND_1'));
print ('B1M', B1M);
print ('B1A', B1A);

var B2M = ee.Image(l8RawSelect.get('RADIANCE_MULT_BAND_2'));
var B2A = ee.Image(l8RawSelect.get('RADIANCE_ADD_BAND_2'));
print ('B2M', B2M);
print ('B2A', B2A);

var B3M = ee.Image(l8RawSelect.get('RADIANCE_MULT_BAND_3'));
var B3A = ee.Image(l8RawSelect.get('RADIANCE_ADD_BAND_3'));
print ('B3M', B3M);
print ('B3A', B3A);

var B4M = ee.Image(l8RawSelect.get('RADIANCE_MULT_BAND_4'));
var B4A = ee.Image(l8RawSelect.get('RADIANCE_ADD_BAND_4'));
print ('B4M', B4M);
print ('B4A', B4A);

var B5M = ee.Image(l8RawSelect.get('RADIANCE_MULT_BAND_5'));
var B5A = ee.Image(l8RawSelect.get('RADIANCE_ADD_BAND_5'));
print ('B5M', B5M);
print ('B5A', B5A);

var B6M = ee.Image(l8RawSelect.get('RADIANCE_MULT_BAND_6'));
var B6A = ee.Image(l8RawSelect.get('RADIANCE_ADD_BAND_6'));
print ('B6M', B6M);
print ('B6A', B6A);

var B7M = ee.Image(l8RawSelect.get('RADIANCE_MULT_BAND_7'));
var B7A = ee.Image(l8RawSelect.get('RADIANCE_ADD_BAND_7'));
print ('B7M', B7M);
print ('B7A', B7A);

var BM = [B1M, B2M, B3M, B4M, B5M, B6M, B7M];
var BA = [B1A, B2A, B3A, B4A, B5A, B6A, B7A];
print (BM);
print (BA);

//create an ee.image to save new radiance bands to
var toaRad = ee.Image(0);

//calculate TOA radiance by looping through bands and applying gain and bias
for (var i=0;i<bands.length; i++){
  
  var band = l8RawSelect.select(bands[i]);
  var mult = ee.Number(BM[i]);
  var add = ee.Number(BA[i]);
  var Rad = (band.multiply(ee.Image(mult)).add(ee.Image(add)));
  
  toaRad = ee.Image(toaRad.addBands(Rad.rename(bands[i])));
}  

//check result
print (toaRad);

//display result
var trueColor2 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 100};
Map.addLayer(toaRad, trueColor2,'Landsat Radiance True Color');  
  
