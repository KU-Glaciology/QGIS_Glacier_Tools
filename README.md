# QGIS_Glacier_Tools
Scripts for use in QGIS for various glaciology related calculations

## Requirements
These scripts are meant to be put into the scripts folder for a local installation of QGIS.  This will make them show up in the QGIS toolbox grouped under the user scripts folder in a Glacier Tools toolbox.

To run you will need to have (rasterio)[https://rasterio.readthedocs.io/en/latest/], (numpy)[https://numpy.org/] and (scipy)[https://www.scipy.org/], I believe installing rasterio will install the others.

QGIS 3.16 includes a version of pip3 that you can use to install python dependencies for use by the QGIS processing toolbox.  See this (issue)[https://github.com/qgis/QGIS-Mac-Packager/issues/33] for some background.  

I have tested on MacOS Mojave using QGIS 3.16.
