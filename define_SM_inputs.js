// Script written by Ryan L. Crumley, July 2019
// Edited by DFH to add PRISM and fix NLCD, July 2019
// Edited by Nina to add SWR & LWR, Dec 2019
 
// When using this script, simply copy and paste to your own file directory system.
// Experiment and make changes to your own version of this script.
// Any changes that you make here will be saved here permanently. Only make changes that work!
// If you make changes, comment them out in a descriptive manner, and leave your name and date in the comments like;
// RLC add, 2019-07-15 Blah Blah Blah


////////////////////////////////////////////////////////////////////////
// This script will create all the required inputs for SnowModel, in geotiff format
// It can be used in conjunction with the Matlab script from D.Hill, July 2019

// OUTPUTS of this script: 
//      1) NLCD of the user-defined region, in geotiff
//      2) DEM of the user-defined region, in geotiff
//      3) Reanlaysis inputs for creating the MicroMet file, in geotiff format
//      4) PRISM climatologies of temp and precip (1 file per month)

// First, define a color scheme for the map visualization below.
var grnbrn = ['#543005','#8c510a','#bf812d','#dfc27d','#f6e8c3','#f5f5f5','#c7eae5','#80cdc1','#35978f','#01665e','#003c30'];
var visparams = {min:0,max:4000,palette:grnbrn};
var visparams_aster = {min:0,max:4000,bands:['elevation'],palette:grnbrn};

//////////////////////////////////////////////////////////////////////
////////////////   Variables requiring input   ///////////////////////
//////////////////////////////////////////////////////////////////////

// Create a domain name to attach to your output. Optional.
var domain_name = 'WY'

// These are the min and max corners of your domain in Lat, Long
// Western Wyoming
// Input the minimum lat, lower left corner
var minLat = 42.363116
// Input the minimum long, lower left corner
var minLong = -111.155208
// Input the max lat, upper right corner
var maxLat = 44.582480
// Input the max Long, upper right corner
var maxLong = -109.477849

// These are the min and max corners of your reanalysis in Lat, Long (create a slightly larger box)
// Input the minimum lat, lower left corner
var minLat2 = (minLat - 0.25);
// print(minLat2);
// Input the minimum long, lower left corner
var minLong2 = (minLong - 0.5);
// Input the max lat, upper right corner
var maxLat2 = (maxLat + 0.25);
// Input the max Long, upper right corner
var maxLong2 = (maxLong + 0.5);

// This resolution for final output of NLCD and DEM only
// This is in meters
var my_resolution = 100

// Define the final output projection using EPSG codes
// WGS UTM Zone 12 Code for Idaho/Wyoming = 32612
// WGS UTM Zone 11 Code for Nevada        = 32611
// WGS UTM Zone 10 Code for West Coast    = 32610
// WGS 84 4326
// WGS UTM 10
// WGS Alaska Albers = 3338
var epsg_code = 'EPSG:32612';

// Name the DEM output
var dem_name = 'DEM';
// Name the Land Cover output
var lc_name = 'NLCD2016';

// The Beginning and End Dates you care about //
// This will start on the 'begin' date at 0:00 and the last iteration will be 
// on the day before the 'end' date below. Look at the printed variable 'tair from CFSv2' 
// in the console to double check.
var begin = '2014-09-01';
var end = '2019-09-01';

//////////////////////////////////////////////////////////////////
/////////////////      DOMAIN     ////////////////////////////////
//////////////////////////////////////////////////////////////////

// Define the desired rectangular domain
// NOTE: The projection is not reset until the exporting process which also
// allows for it to be visualized as a layer.
var my_domain = ee.Geometry.Rectangle({
  coords:[minLong,minLat,maxLong,maxLat],
  proj: 'EPSG:4326',
  geodesic:true,
});
//,'EPSG:32612',false);

// This adds the domain you care about to the visualization
Map.addLayer(my_domain,visparams,'My Domain');
print (my_domain);

// This adds the extent of the reanalysis product to the visualization.
var my_domain2 = ee.Geometry.Rectangle([minLong2,minLat2,maxLong2,maxLat2]);//,'EPSG:32612',false);
Map.addLayer(my_domain2,visparams,'My Reanalysis Domain');
print (my_domain2);

// Check the domain area in meters squared. Uncomment to check.
//var my_domain_area = my_domain.area();
//print(my_domain_area);


////////////////   Datasets of Interest  //////////////////////
////////    Digital Elevation Models and Land Cover   /////////
///////////////////////////////////////////////////////////////

// NOTE: several choices below for DEM. Uncomment your preferred option

////////   Import 30m SRTM Data   ///////////////////
// NOTE: This only covers through 60 degrees latitude. See visualization layers.
//var SRTM30 = ee.Image('USGS/SRTMGL1_003');
// Find out some info about this image (hint: look in the console)
//var bands30 = SRTM30.bandNames();
//var info30 = SRTM30.getInfo();
//print(bands30,'Band Names');
//print(info30,'Band Info');
//Map.addLayer(SRTM30,visparams,'SRTM30');

////////  Import 100m ASTER data //////////////
// NOTE: this works above 60 deg lat; better for Alaska...
//var ASTER = ee.Image('NASA/ASTER_GED/AG100_003');
// Find out some info about this image (hint: look in the console)
//var bands100 = ASTER.bandNames();
//var info100 = ASTER.getInfo();
//print(bands100,'Band Names');
//print(info100,'Band Info');
//Map.addLayer(ASTER,visparams_aster,'ASTER');

/////////  Import 90m SRTM Data   ////////////////////
// NOTE: This only covers through 60 degrees latitude. See visualization layers.
var SRTM90 = ee.Image('CGIAR/SRTM90_V4');
var bands90 = SRTM90.bandNames();
var info90 = SRTM90.getInfo();
print(bands90,'Band Names');
print(info90,'Band Info');
//Map.addLayer(SRTM90,visparams,'SRTM90');

////////   Import NLCD Dataset   ////////////////////
var NLCD = ee.ImageCollection('USGS/NLCD');
//Next: the NLCD has numerous images for different years. I want to use
//the most current (2016), so I filter by time to isolate 2016 slice.
var landcover = NLCD.select('landcover');
var landcoverfiltered=landcover.filterDate('2015-01-01','2018-01-01');
var landcoverVis = {
  min: 0.0,
  max: 95.0,
  palette: [
    '466b9f', 'd1def8', 'dec5c5', 'd99282', 'eb0000', 'ab0000', 'b3ac9f',
    '68ab5f', '1c5f2c', 'b5c58f', 'af963c', 'ccb879', 'dfdfc2', 'd1d182',
    'a3cc51', '82ba9e', 'dcd939', 'ab6c28', 'b8d9eb', '6c9fb8'
  ],
};
// Next, the following is the way I have been able to convert the image collection 
// (only one image at this point) to a single image. Must be a better way.
var lcsingle=landcoverfiltered.median();
Map.addLayer(lcsingle, landcoverVis, 'Landcover');

////////////////   Datasets of Interest  //////////////////////
////////                PRISM DATA                    /////////
///////////////////////////////////////////////////////////////

//NOTE: these are 2.5 arc min. Roughly ~4 km. Seems to be the best option available...
//NOTE: you can pick any 30 year period you want. I chose 1985-2015.

/////////  Import PRISM Climatologies   ////////////////////
// Precip first...
var prism = ee.ImageCollection('OREGONSTATE/PRISM/AN81m');
var precipitation = prism.select('ppt');
var janppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(1,1,'month'));
var febppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(2,2,'month'));
var marppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(3,3,'month'));
var aprppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(4,4,'month'));
var mayppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(5,5,'month'));
var junppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(6,6,'month'));
var julppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(7,7,'month'));
var augppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(8,8,'month'));
var sepppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(9,9,'month'));
var octppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(10,10,'month'));
var novppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(11,11,'month'));
var decppt = precipitation.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(12,12,'month'));

// reduce image collection with mean()
var janmean = janppt.mean();
var febmean = febppt.mean();
var marmean = marppt.mean();
var aprmean = aprppt.mean();
var maymean = mayppt.mean();
var junmean = junppt.mean();
var julmean = julppt.mean();
var augmean = augppt.mean();
var sepmean = sepppt.mean();
var octmean = octppt.mean();
var novmean = novppt.mean();
var decmean = decppt.mean();

var precipitationVis = {
  min: 0.0,
  max: 300.0,
  palette: ['red', 'yellow', 'green', 'cyan', 'purple'],
};
Map.setCenter(-100.55, 40.71, 4);
Map.addLayer(janmean, precipitationVis, 'Jan Precipitation');
Map.addLayer(febmean, precipitationVis, 'Feb Precipitation');
Map.addLayer(marmean, precipitationVis, 'Mar Precipitation');
Map.addLayer(aprmean, precipitationVis, 'Apr Precipitation');
Map.addLayer(maymean, precipitationVis, 'May Precipitation');
Map.addLayer(junmean, precipitationVis, 'Jun Precipitation');
Map.addLayer(julmean, precipitationVis, 'Jul Precipitation');
Map.addLayer(augmean, precipitationVis, 'Aug Precipitation');
Map.addLayer(sepmean, precipitationVis, 'Sep Precipitation');
Map.addLayer(octmean, precipitationVis, 'Oct Precipitation');
Map.addLayer(novmean, precipitationVis, 'Nov Precipitation');
Map.addLayer(decmean, precipitationVis, 'Dec Precipitation');

// tmean next
var tmean = prism.select('tmean');
var jantmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(1,1,'month'));
var febtmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(2,2,'month'));
var martmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(3,3,'month'));
var aprtmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(4,4,'month'));
var maytmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(5,5,'month'));
var juntmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(6,6,'month'));
var jultmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(7,7,'month'));
var augtmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(8,8,'month'));
var septmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(9,9,'month'));
var octtmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(10,10,'month'));
var novtmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(11,11,'month'));
var dectmean = tmean.filter(ee.Filter.calendarRange(1985,2015,'year'))
.filter(ee.Filter.calendarRange(12,12,'month'));

// reduce image collection with mean()
var jantmean = jantmean.mean();
var febtmean = febtmean.mean();
var martmean = martmean.mean();
var aprtmean = aprtmean.mean();
var maytmean = maytmean.mean();
var juntmean = juntmean.mean();
var jultmean = jultmean.mean();
var augtmean = augtmean.mean();
var septmean = septmean.mean();
var octtmean = octtmean.mean();
var novtmean = novtmean.mean();
var dectmean = dectmean.mean();

var tmeanVis = {
  min: -30.0,
  max: 30.0,
  palette: ['red', 'yellow', 'green', 'cyan', 'purple'],
};
Map.setCenter(-100.55, 40.71, 4);
Map.addLayer(jantmean, tmeanVis, 'Jan Tmean');
Map.addLayer(febtmean, tmeanVis, 'Feb Tmean');
Map.addLayer(martmean, tmeanVis, 'Mar Tmean');
Map.addLayer(aprtmean, tmeanVis, 'Apr Tmean');
Map.addLayer(maytmean, tmeanVis, 'May Tmean');
Map.addLayer(juntmean, tmeanVis, 'Jun Tmean');
Map.addLayer(jultmean, tmeanVis, 'Jul Tmean');
Map.addLayer(augtmean, tmeanVis, 'Aug Tmean');
Map.addLayer(septmean, tmeanVis, 'Sep Tmean');
Map.addLayer(octtmean, tmeanVis, 'Oct Tmean');
Map.addLayer(novtmean, tmeanVis, 'Nov Tmean');
Map.addLayer(dectmean, tmeanVis, 'Dec Tmean');


////////////////   Datasets of Interest  //////////////////////
////////                Reanalysis DATA               /////////
///////////////////////////////////////////////////////////////
var cfsv2 = ee.ImageCollection('NOAA/CFSV2/FOR6H')
                  .filter(ee.Filter.date(begin,end));
var tair = cfsv2.select('Temperature_height_above_ground').toBands();
var elev = cfsv2.select('Geopotential_height_surface').toBands();
var uwind = cfsv2.select('u-component_of_wind_height_above_ground').toBands();
var vwind = cfsv2.select('v-component_of_wind_height_above_ground').toBands();
var surfpres = cfsv2.select('Pressure_surface').toBands();
var spechum = cfsv2.select('Specific_humidity_height_above_ground').toBands();
var prec = cfsv2.select('Precipitation_rate_surface_6_Hour_Average').toBands();
var lwr = cfsv2.select('Downward_Long-Wave_Radp_Flux_surface_6_Hour_Average').toBands();
var swr = cfsv2.select('Downward_Short-Wave_Radiation_Flux_surface_6_Hour_Average').toBands();
// To check the time iterations, look at the printed variable in the console
print(tair, 'tair from CFSv2');

//////////////////////////////////////////////////////////////
///////   EXPORT, RESCALE, REPROJECT, CLIP  //////////////////
//////////////////////////////////////////////////////////////

// Export the SRTM DEM to Geotiff 
//Export.image.toDrive({
//  image: SRTM90,
//  description: dem_name+'_'+domain_name,
//  region: my_domain,
//  scale: my_resolution,
//  crs: epsg_code,
//  maxPixels: 1e12
//});

// Export the DEM to Geotiff 
Export.image.toDrive({
  image: SRTM90,
  description: dem_name+'_'+domain_name,
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
  maxPixels: 1e12
});

// Export the NLCD to Geotiff
Export.image.toDrive({
  image: lcsingle,
  description: lc_name+'_'+domain_name,
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the CFSv2 Temp to Geotiff 
Export.image.toDrive({
  image: tair,
  description: 'cfsv2_'+begin+end+'_tair',
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});


// Export the CFSv2 Temp to Geotiff 
Export.image.toDrive({
  image: elev,
  description: 'cfsv2_'+begin+end+'_elev',
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});

// Export the CFSv2 Prec to Geotiff 
Export.image.toDrive({
  image: prec,
  description: 'cfsv2_'+begin+end+'_prec' ,
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});

// Export the CFSv2 Uwind to Geotiff 
Export.image.toDrive({
  image: uwind,
  description: 'cfsv2_'+begin+end+'_uwind' ,
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});


// Export the CFSv2 Vwind to Geotiff 
Export.image.toDrive({
  image: vwind,
  description: 'cfsv2_'+begin+end+'_vwind' ,
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});

// Export the CFSv2 Surfpres to Geotiff 
Export.image.toDrive({
  image: surfpres,
  description: 'cfsv2_'+begin+end+'_surfpres' ,
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});

// Export the CFSv2 RedHum to Geotiff 
Export.image.toDrive({
  image: spechum,
  description: 'cfsv2_'+begin+end+'_spechum' ,
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});

// Export the CFSv2 LWR to Geotiff 
Export.image.toDrive({
  image: lwr,
  description: 'cfsv2_'+begin+end+'_lwr' ,
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});

// Export the CFSv2 SWR to Geotiff 
Export.image.toDrive({
  image: swr,
  description: 'cfsv2_'+begin+end+'_swr' ,
  region: my_domain2,
  scale: 22200,
  crs: epsg_code,
});

// Precip grids
// Export the Jan ppt to Geotiff
Export.image.toDrive({
  image: janmean,
  description: 'janppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Feb ppt to Geotiff
Export.image.toDrive({
  image: febmean,
  description: 'febppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Mar ppt to Geotiff
Export.image.toDrive({
  image: marmean,
  description: 'marppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Apr ppt to Geotiff
Export.image.toDrive({
  image: aprmean,
  description: 'aprppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the May ppt to Geotiff
Export.image.toDrive({
  image: maymean,
  description: 'mayppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Jun ppt to Geotiff
Export.image.toDrive({
  image: junmean,
  description: 'junppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Jul ppt to Geotiff
Export.image.toDrive({
  image: julmean,
  description: 'julppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Aug ppt to Geotiff
Export.image.toDrive({
  image: augmean,
  description: 'augppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Sep ppt to Geotiff
Export.image.toDrive({
  image: sepmean,
  description: 'sepppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Oct ppt to Geotiff
Export.image.toDrive({
  image: octmean,
  description: 'octppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Nov ppt to Geotiff
Export.image.toDrive({
  image: novmean,
  description: 'novppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Dec ppt to Geotiff
Export.image.toDrive({
  image: decmean,
  description: 'decppt',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Tmean grids
// Export the Jan tmean to Geotiff
Export.image.toDrive({
  image: jantmean,
  description: 'jantmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Feb tmean to Geotiff
Export.image.toDrive({
  image: febtmean,
  description: 'febtmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Mar tmean to Geotiff
Export.image.toDrive({
  image: martmean,
  description: 'martmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Apr tmean to Geotiff
Export.image.toDrive({
  image: aprtmean,
  description: 'aprtmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the May tmean to Geotiff
Export.image.toDrive({
  image: maytmean,
  description: 'maytmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Jun tmean to Geotiff
Export.image.toDrive({
  image: juntmean,
  description: 'juntmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Jul tmean to Geotiff
Export.image.toDrive({
  image: jultmean,
  description: 'jultmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Aug tmean to Geotiff
Export.image.toDrive({
  image: augtmean,
  description: 'augtmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Sep tmean to Geotiff
Export.image.toDrive({
  image: septmean,
  description: 'septmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Oct tmean to Geotiff
Export.image.toDrive({
  image: octtmean,
  description: 'octtmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Nov tmean to Geotiff
Export.image.toDrive({
  image: novtmean,
  description: 'novtmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});

// Export the Dec tmean to Geotiff
Export.image.toDrive({
  image: dectmean,
  description: 'dectmean',
  region: my_domain,
  scale: my_resolution,
  crs: epsg_code,
});


///////////////////////////////////////////////////////////////
/////////////////  EXTRA STUFF ////////////////////////////////
///////////////////////////////////////////////////////////////

/*
var dataset = ee.Image('JAXA/ALOS/AW3D30_V1_1');
var elevation = dataset.select('AVE');
var elevationVis = {
  min: 0.0,
  max: 4000.0,
  palette: ['0000ff', '00ffff', 'ffff00', 'ff0000', 'ffffff'],
};
Map.setCenter(136.85, 37.37, 4);
Map.addLayer(elevation, elevationVis, 'Elevation');
*/

///////////////////////////////////////////////////////////////
////////////     Import HUC Watersheds     ////////////////////
///////////////////////////////////////////////////////////////

// If you want to use the HUC watersheds for some reason, uncomment these lines.

var HUC = ee.FeatureCollection('USGS/WBD/2017/HUC08');
var styleParams = {
  fillColor: 'ece7f2',
  color: '000000',
  width: 1.0,
};
var wsheds = HUC.style(styleParams);
Map.addLayer(wsheds, {}, 'USGS/WBD/2017/HUC08');



////////////////////////////////////////////////////////////////
/////////  Save Info from domains in the script    /////////////
////////////////////////////////////////////////////////////////

/*
// Thompson Pass Domain
// Input the minimum lat, lower left corner
var minLat = 60.9651385
// Input the minimum long, lower left corner
var minLong = -146.4828057
// Input the max lat, upper right corner
var maxLat = 61.538588
// Input the max Long, upper right corner
var maxLong = -144.879882
//GOA Domain
// Input the minimum lat, lower left corner
var minLat = 56.2819
// Input the minimum long, lower left corner
var minLong = -156.8955
// Input the max lat, upper right corner
var maxLat = 60.8622
// Input the max Long, upper right corner
var maxLong = -122.7722
*/

//Central OR Domain
// Input the minimum lat, lower left corner
//var minLat = 42.045789
// Input the minimum long, lower left corner
//var minLong = -123.476288
// Input the max lat, upper right corner
//var maxLat = 45.702675
// Input the max Long, upper right corner
//var maxLong = -121.231792


//////////////////////////////////////////////////////////
/////////////   ASPECT/SLOPE/HILLSHADE  //////////////////
///////////////////////////////////////////////////////////////
/*
// This pre-cooked GEE function calculates slope in degrees (0-90) from the DEM layer.
// It uses a 4 connected neighbors approach and edge pixels will have missing data.
var slope = ee.Terrain.slope(DEM);
Map.addLayer(slope,{},'Slope');
// This pre-cooked GEE function calculates aspect in degrees (0-365) from the DEM layer.
// It uses a 4 connected neighbors approach and edge pixels will have missing data.
var aspect = ee.Terrain.aspect(DEM);
Map.addLayer(aspect,{},'Aspect');
// This pre-cooked GEE function creates a hillshade layer from the DEM.
// This can come in handy when trying to recognize local geography and physical features.
var hillshade = ee.Terrain.hillshade(DEM);
Map.addLayer(hillshade,visHillshade,'Hillshade');
// This pre-cooked function allows all of the terrain products to be visualized in a single,
// multi-band image. I like the red tinted visualization the best, with hillshade, slope, and
// elevation data as the RGB layers. 
var all = ee.Terrain.products(DEM);
Map.addLayer(all,visProducts,'All Terrain Products');
*/
