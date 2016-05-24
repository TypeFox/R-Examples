# file example.KML.GeotopRasterBrick.R 
#
# This file contains a script which plots the measurement sites on a GoogleMap support and uses the function 'brickFromOutputSoil3DTensor'
# WARNING: the simulation template contanis only few maps utilized in this script in order to save memory disk
# The data and te results contained in this script are for educational use only and may not be realistic.
#
#
# author: Emanuele Cordano on 03-07-2013
#
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

library(geotopbricks)

wpath <- "http://www.rendena100.eu/public/geotopbricks/simulations/idroclim_test1"

x <- "SoilLiqContentTensorFile"
when <- as.POSIXct("2002-03-22 UTC",tz="Etc/GMT+1")
b <- brickFromOutputSoil3DTensor(x,when=when,wpath=wpath,tz="A",use.read.raster.from.url=TRUE)

# set colors 
N <- 1000
start_col <- 4/12
end_col <- 9/12
col <- rainbow(start=start_col,end=end_col,n=N,alpha=0.6)

# FILE kmz

wpath_kmz <- '.'

filename_kmz <- paste(x,as.character(when),sep="_")
filename_kmz <- str_replace(filename_kmz," ","_")
filename_kmz <- paste(wpath_kmz,filename_kmz,sep="/")
filename_kmz <- paste(filename_kmz,".kmz",sep="")

KML(geotopbrick(b),filename=filename_kmz,col=col,overwrite=TRUE)








