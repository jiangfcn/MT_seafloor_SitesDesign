# MTsitesDesign
 MATLAB package for plotting a seafloor MT site profile with water depths.

# Notes:
 * GMT (the Generic Mapping Tools) and Google Earth should be installed before using this package for MT site design.
 * A DEM file that GMT functions can read is needed.

# Guidelines/Workflow
 1. Copy all the scripts into a folder and add it to your MATLAB PATH.
 2. Generate points with equal spacing on a profile by using PointsOnGeolines.m, open up this script, and edit the parameters 'n', 'filename', 'ID', 'P1', and 'P2'. What these parameters indicate have been annotated on each line.
 3. Input the locations file to Google Earth, and edit points by hand. Then, output a kml file which contains the final MT site locations you want.
 4. Copy InputParameters.m to your work folder and open it, editing the structure array u.
 5. Run TopoTrackProfile, and a MT site profile with topographic variations should be generated automatically.


 # Copyright
  * Copyright 2022-2023 Feng Jiang South China Sea Institute of Oceanology.
  * Functions utm2geo and geo2utm are included in this package, shared or edited by Dennis Rippe, 2010. Please contact the author if there are copyright problems.

 # License
  * This package is free to use under the GNU General Public License terms.
  * GMT functions are called in part of scripts. License corresponding to GMT, if any, should be complied with. http://gmt.soest.hawaii.edu/projects/gmt
  
