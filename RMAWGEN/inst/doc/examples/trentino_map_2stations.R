# file trentino_map_2stations.R
#
# This file contains a script which plots the measurement sites on a GoogleMap support
#
#
# author: Emanuele Cordano on 09-08-2012

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





# This script plots the measurement sites on a GoogleMap support 

rm(list=ls())
library(RMAWGEN)
library(RgoogleMaps)
data(trentino)

TEMP <- TEMPERATURE_MAX[,-(1:3)]

# set the 2 stations
station <- c("T0090","T0083")



coord <- STATION_LATLON[STATION_NAMES  %in% names(TEMP),c(1,2)]
rownames(coord) <- STATION_NAMES[STATION_NAMES  %in% names(TEMP)]

iday <- extractdays(when="1985-7-1",origin="1958-1-1")
TEMPv <- TEMP[iday,]
col <- array(NA,length(TEMPv))

# set colors 
start <- 4/6 
end=0.99999
TEMP_max <- max(TEMP,na.rm=TRUE)
TEMP_min <- min(TEMP,na.rm=TRUE)

col_se <- (range(TEMPv,na.rm=TRUE)-TEMP_min)/(TEMP_max-TEMP_min)*(end-start)+start

names(col) <- names(TEMPv)
col[names(sort(TEMPv))] <- rainbow(start=min(col_se),end=max(col_se),n=length(sort(TEMPv))) # color scale from "blue" to "red"

# set bounding box 
lonR <- range(coord[,1])
latR <- range(coord[,2])

col2 <- col[station]
coord2 <- coord[col %in% col2,]
TrentinoMap <- GetMap.bbox(lonR=lonR,latR=latR)
map2 <- PlotOnStaticMap(TrentinoMap,lat=coord2[station,2],lon=coord2[station,1],FUN=points,cex=1.5,pch=20,col=col2,add=FALSE)

#TextOnStaticMap(map2,lat=coord2[station,2],lon=coord2[station,1],labels=station, offset=1,add=TRUE)
# end of the script

