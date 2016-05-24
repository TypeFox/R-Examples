# file example.KML.GeotopRasterBrick.R 
#
# This file contains a script which plots the measurement sites on a GoogleMap support
#
#
# author: Emanuele Cordano on 08-11-2012

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
library(rgdal)
library(raster)
library(zoo)
library(geotopbricks)

## working path for 3D distributed raster maps
## The study case of Ton-Toss (Val di Non, Trentino, Italy)

wpath <- 'http://meteogis.fmach.it/idroclima/ton-toss'
# WARNING: In order to save disk space, some files of this simulation (unusuful for the example) were removed !!!!
# keyword for water content 3D+time output raster maps
watercontent_prefix <- get.geotop.inpts.keyword.value("SoilLiqContentTensorFile",wpath=wpath) #"thetaliq"
#  crs projection
crs  <-"+proj=utm +zone=32 +ellps=WGS84 +datum=WGS84 +units=m +no_defs "

#  vector with vertical layer thckness

layers <- get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath)
names(layers) <- paste("L",1:length(layers))

# set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")

start <-  get.geotop.inpts.keyword.value("InitDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")
end <- get.geotop.inpts.keyword.value("EndDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")

# set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")


# In this examples maps are provided with daily frequency!!
time <- seq(from=start,to=end,by="days")

#
# files extracts the filename of the maps of the first 4 layers!!
#
files <- pointer.to.maps.xyz.time(wpath=wpath,map.prefix=watercontent_prefix,zoo.index=time,nlayers=4)

# import maps

start_m <- as.POSIXlt("2012-04-15 00:00",tz="A")
end_m <- as.POSIXlt("2012-04-20 00:00",tz="A")

gt_wtc1 <- geotopbrick(x=files,layer=1,timerange=c(start_m,end_m),crs=crs)
gt_wtc2 <- geotopbrick(x=files,layer=2,timerange=c(start_m,end_m),crs=crs)
gt_wtc3 <- geotopbrick(x=files,layer=3,timerange=c(start_m,end_m),crs=crs)

# Averaged Soil Water Content in the first 3 soil layers (about 33 centimetes)
# 'RasterBrick objects'
wtc <- (brick(gt_wtc1)*layers[1]+brick(gt_wtc2)*layers[2]+brick(gt_wtc3)*layers[3])/sum(layers[1:3])
# 'GeotopRasterBrick objects
gt_wtc <- ((gt_wtc1)*layers[1]+(gt_wtc2)*layers[2]+(gt_wtc3)*layers[3])/sum(layers[1:3])

# set colors for 'KML' and 'plot'

N <- 10000

start_watercontent <- 4/12
end_watercontent <- 9/12
col_watercontent <- rainbow(start=start_watercontent,end=end_watercontent,n=N,alpha=0.6)
# Creates the KML
KML(gt_wtc,filename="zz_wtc_33cm_ton_toss.kml",overwrite=TRUE,col=col_watercontent)

# raster legend (in
pdf <- "zz_wtc_33cm_ton_toss_legend.pdf"
color.bar.raster(x=gt_wtc,col_watercontent,digits=2,pdf=pdf)