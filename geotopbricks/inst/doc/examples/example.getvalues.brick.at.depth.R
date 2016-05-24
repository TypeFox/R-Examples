# file example.getvalues.brick.at.depth.R
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







library(rgdal)
library(raster)
library(zoo)
library(methods)
library(geotopbricks)
wpath <- 'http://meteogis.fmach.it/idroclima/panola13_run2xC_test3'

prefix <- get.geotop.inpts.keyword.value("SoilLiqWaterPressTensorFile",wpath=wpath)


slope <- get.geotop.inpts.keyword.value("SlopeMapFile",raster=TRUE,wpath=wpath)
bedrock_depth <- get.geotop.inpts.keyword.value("BedrockDepthMapFile",raster=TRUE,wpath=wpath)



layers <- get.geotop.inpts.keyword.value("SoilLayerThicknesses",numeric=TRUE,wpath=wpath)

names(layers) <- paste("L",1:length(layers),sep="")

# set van genuchten parameters to estimate water volume
theta_sat <- get.geotop.inpts.keyword.value("ThetaSat",numeric=TRUE,wpath=wpath)
theta_res <- get.geotop.inpts.keyword.value("ThetaRes",numeric=TRUE,wpath=wpath)
alphaVG <-  get.geotop.inpts.keyword.value("AlphaVanGenuchten",numeric=TRUE,wpath=wpath) # expressed in mm^-1
nVG <-  get.geotop.inpts.keyword.value("NVanGenuchten",numeric=TRUE,wpath=wpath)


# end set van genuchten parameters to estimate water volume


# set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")

start <-  get.geotop.inpts.keyword.value("InitDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")
end <- get.geotop.inpts.keyword.value("EndDateDDMMYYYYhhmm",date=TRUE,wpath=wpath,tz="A")

# end set time during wich GEEOtop simulation provided maps (the script is written for daily frequency")

## The maps are obtanied with daily frequancy
time <- seq(from=start,to=end,by="days")

## Pressure head map filename

psiliq_files <- pointer.to.maps.xyz.time(wpath=wpath,map.prefix=prefix,zoo.index=time,nlayers=length(layers))

## Note that in this similation 'psi' maps are returned with daily frequency!!!!
days <- c(1,2,5,10,15,20) # integer!!!

i <- 5 # ten days since the beginning of the simulation period

psi <- brick(psiliq_files,layer=1:length(layers),rows=days[i])
psi_bedrock <- getvalues.brick.at.depth(psi,depth=bedrock_depth,layers=layers,verify=FALSE)
# Plot the map of the presuure head over the bedrock:
# plot(psi_bedrock$val_z0) # expressed in millimeters !!!! - Do not Run

## Calculation of water-table thickness according to Cordano and Rigon's Equation 10 (http://www.agu.org/pubs/crossref/2008/2006WR005740.shtml)

slope_cos <- cos(slope*pi/180)

## geotop_watertable_thickness: i.e. the slope-normal thickness of water-table film over the bedrock!!!

geotop_watertable_thickness <- (psi_bedrock$val_z0/slope_cos+bedrock_depth-psi_bedrock$z0)
geotop_watertable_thickness[geotop_watertable_thickness<0] <- 0

# plot(geotop_watertable_thickness) # expressed in millimeters !!!! - Do not Run