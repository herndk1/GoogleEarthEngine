var lstoaCollection = ee.ImageCollection("LANDSAT/LC08/C01/T1_RT_TOA"),
    pt = /* color: #d63000 */ee.Geometry.Point([5.02899169921875, 15.392783791106002]);


//identify LS TOA reflectance scene
//filter landsat TOA reflectance image collection
var lstoa = ee.Image(lstoaCollection
  .filterBounds(pt)
  .filterDate('2015-10-20', '2015-10-22')
  .first());
  
//display resulting Landsat ToA reflectance scene
var trueColor1 = {bands: ['B4', 'B3', 'B2'], min: 0, max: 0.3};
Map.addLayer(lstoa, trueColor1,'Landsat TOA reflectance true color');


//define function to apply dark object subtraction method
function dos(img) {
  //define bands to apply dos 
  var bands = ['B2','B3','B4','B5','B6','B7'];
  //create ee.Image to add dos corrected bands to
  var sr = ee.Image(0);
  //create loop to apply DOS to each band
  for (var i=0;i<bands.length;i++){
    //identify the band you will be applying the dos to
    var band = ee.Image(img.select(bands[i]).rename('value'));
   
    //identify the value for the darkest 1% of pixels or the minimum pixel (choose 1 of the lines below)
    var darkObject = band.reduceRegion({reducer:ee.Reducer.percentile([1]),geometry:img.geometry(),scale:150,bestEffort:true});
    //var darkObject = band.reduceRegion({reducer:ee.Reducer.min(),geometry:img.geometry(),scale:150,bestEffort:true});
    
    //get the value of the darkest 1%
    var doValue = ee.Number(darkObject);
    print (doValue);
    //make DO value an ee.Image
    var doImg = ee.Image((ee.Number(darkObject.get('value'))));
    //subtract the DO value from the entire scene
    var srband = band.subtract(doImg);
    //add the corrected bands to the ee.Image 
    sr = sr.addBands(srband.rename(bands[i]));
  }
  return sr;
}

//apply dos to selected landsat scene
var lsDOS = dos(lstoa);

//display atmospherically corrected Landsat 8 sr image
Map.addLayer(lsDOS, trueColor1, 'L8 DOS correction surface reflectance')
