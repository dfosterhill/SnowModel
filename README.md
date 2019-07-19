# Overview

This repository stores for scripts related to pre- / post-processing files for [Glen Liston's](https://www.cira.colostate.edu/staff/liston-glen/) MicroMet/SnowModel/HydroFlow models. See manuscript [here] for more details on the modeling suite. 

## Matlab

1. GEE_to_snowmodel.m - this script will take output from Google Earth Engine and it will convert it to file formats needed for SnowModel. Please see the comments at the top of that file. 
2. arcgridwrite.m - this script is called by GEE_to_snowmodel.m for the purposes of writing out ESRI ASCII grids.

## Google Earth Engine

1. define_sm_inputs.js

You will need to be registered with a Google Earth Engine account to run these scripts. If you don't already have an account you can sign up for one [here](https://signup.earthengine.google.com/#!/).

1. In GEE, in the left toolbar, select 'New' - 'File' from the red dropdown box. Enter a path, a filename, and a description.
2. In the code editor (upper center panel of GEE) paste the contents from the define_sm_inputs.js file provided here in Github. Click on 'Save'
3. Click on 'Run'. At this point you can explore the three tabs in the right toolbar (Inspecctor / Console / Tasks). Additionally, you will see data layers begin to appear in the bottom panel of GEE.
4. Choose 'Tasks'. You will see many image names with a blue 'run' button to the right of each. For each of these, you must click 'run'. This will bring up a dialog box asking if you want to initiate export. You can hit return and it will do so. Repeat this for all images.
5. Go to wherever you directed output (you have that option in the dialog box). Most commonly, this will be your Google Drive. You will find all exported geotiff images there.

## contact 

[David Hill](mailto:dfh@oregonstate.edu)
