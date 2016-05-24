# file example.meteostation.R
#
# This file contains ...
#
# author: Emanuele Cordano on 03-07-2013

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################


rm(list=ls())
library(stringr)
library(sp)
library(rgdal)
library(geotopbricks)
library(plotKML)
wpath <- "http://www.rendena100.eu/public/geotopbricks/simulations/idroclim_test1"
wpath_kml <- '.'
filename_kml <- 'meteopoint.kml'
filename_kml <- paste(wpath_kml,filename_kml,sep="/")

# Example from geotop.ints file 

#NumberOfMeteoStations=20
#MeteoStationCoordinateX=650848.94,667418.75,644512.25,673925.5,658436.12,656731.75,663628.5,691493,646362.06,657311.62,656926.31,664587.69,691634.94,621316.62,653434.94,663596.25,633824.38,669940.5,700104,650793.5
#MeteoStationCoordinateY=5083241,5117975.5,5097349,5110493,5125738,5082504,5116620.5,5102047,5085824,5072279.5,5136199.5,5098670.5,5103481,5078208.5,5094289,5144884,5083901,5096203,5130422.5,5137986
#MeteoStationElevation=957,696,492,983,324,170,203,420,84,172,656,185,412,385,552,918,705,706,1001,773
#MeteoStationWindVelocitySensorHeight=10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10
#MeteoStationTemperatureSensorHeight=2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2,2
#MeteoStationLatitude=45.886196,46.194763,46.014465,46.125839,46.266727,45.87825,46.183498,46.045193,45.9104,45.786148,46.361172,46.021824,46.058048,45.84676,45.985001,46.437698,45.895657,45.998322,46.297848,46.378632
#MeteoStationLongitude=10.944168,11.169777,10.866802,11.251293,11.056045,11.0197,11.12022,11.47508,10.88716,11.02383,11.039959,11.126421,11.47749,10.562442,10.981014,11.129643,10.725044,11.194649,11.59826,10.960859

MeteoStationLatitude <- get.geotop.inpts.keyword.value("MeteoStationLatitude",numeric=TRUE,wpath=wpath)
MeteoStationLongitude <- get.geotop.inpts.keyword.value("MeteoStationLongitude",numeric=TRUE,wpath=wpath)


meteopoints <- data.frame(lat=MeteoStationLatitude,lon=MeteoStationLongitude,meteo_geotop_id=sprintf("id%04d",1:length(MeteoStationLatitude)),stringsAsFactors=FALSE)

coordinates(meteopoints) <- ~lon+lat
projection(meteopoints) <- CRS(as.character("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"))

kml(meteopoints,file=filename_kml,points_names=meteopoints@data$meteo_geotop_id,alpha = 0.8)




