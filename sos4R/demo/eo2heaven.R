# Copyright (C) 2010 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: 	Daniel Nuest (daniel.nuest@uni-muenster.de)
#			Edzer Pebesma (edzer.pebesma@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")

###############################################################################
# EO2Heaven SOS @ TU Dresden
eo2h <- SOS(url = "http://141.30.100.135:8080/eo2heavenSOS/sos")
eo2h
sosContents(eo2h)

###############################################################################
# O2HEAVEN (Earth Observation and Environmental Modelling for the Mitigation of
# Health Risks) is a research project co-funded by the European Commission as
# part of the 7th Framework Programme (FP7) Environmental theme.
#
# It started on 1st February 2010. EO2HEAVEN contributes to a better
# understanding of the complex relationships between environmental changes and
# their impact on human health.
###############################################################################

###############################################################################
# NOX data for Durban, South Africa
###############################################################################
offering_nox <- sosOfferings(eo2h)[["NOX"]]
sosObservedProperties(offering_nox)
sosFeaturesOfInterest(offering_nox)

#observations_nox <- getObservation(sos = eo2h, offering = offering_nox,
#		eventTime = sosCreateTime(sos = eo2h, time = "2011-01-01::2011-01-03"))
# 29 warnings, because of missing parsers:
#warnings()

#eo2h_converters <- SosDataFieldConvertingFunctions(
#		"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString,
#		"http://www.opengis.net/def/property/OGC/0/SamplingTime" = sosConvertTime)
# FIELDS ADDED IN VERSION 0.2-9
#eo2h <- SOS(url = sosUrl(eo2h))
#		dataFieldConverters = eo2h_converters)

observations_nox <- getObservation(sos = eo2h, offering = offering_nox,
		inspect = TRUE,
		eventTime = sosCreateTime(sos = eo2h, time = "2011-01-01::2011-01-03"))
#[sos4R] Received response (size: 56288 bytes), parsing ...
#[sos4R] Finished getObservation to http://141.30.100.135:8080/eo2heavenSOS/sos 
#--> received 7 observation(s) having 315 result values [ 48, 48, 46, 41, 48, 48, 36 ]. 
str(observations_nox, max.level = 5)
sosBoundedBy(observations_nox)
result_nox_bbox <- sosBoundedBy(observations_nox, bbox = TRUE) # sp ready!
plot(result_nox_bbox)

# result extraction
result_nox <- sosResult(observations_nox)
summary(result_nox)

observations_nox[1:2] # subsetting
result_nox_1_2 <- sosResult(observations_nox[1:2])
summary(result_nox_1_2)
# one feature of interest per result observatoin/per member

result_nox_1 <- sosResult(observations_nox[[1]])
summary(result_nox_1)
names(result_nox_1)
str(observations_nox[[1]])


###############################################################################
# create plot
plotStationNOX <- function(observationId, add = FALSE, colour, ...) {
	.x <- sosResult(observations_nox[[observationId]])[["SamplingTime"]]
	.y <- sosResult(observations_nox[[observationId]])[["NOX"]]
	if(add) {
		lines(x = .x, y = .y, col = colour[[observationId]], ...)
	}
	else {
		plot(x = .x,
				y = .y,
				type = "l", # "b",
				col = colour[[observationId]],
				...,
		)
	}
}

library("RColorBrewer") # display.brewer.all()
stations_colours <- brewer.pal(8, "Dark2")

plotStationNOX(observationId = 1, colour = stations_colours, ylim = c(0, 120),
		xlab = NA, ylab = NA)
plotStationNOX(observationId = 2, colour = stations_colours, add = TRUE)
plotStationNOX(observationId = 3, colour = stations_colours, add = TRUE)
plotStationNOX(observationId = 4, colour = stations_colours, add = TRUE)
plotStationNOX(observationId = 5, colour = stations_colours, add = TRUE)

stations_features <- unlist(sosFeatureIds(observations_nox[1:5]))
stations_features_short <- lapply(stations_features, FUN = substring, first = 63)
title(main = paste("Features:", toString(stations_features_short)),
		xlab = "sampling time", ylab = "NOX")
legend("topleft",  
		legend = stations_features_short,  
		col=stations_colours,  
		lty=1, lwd=1,
		title="NOX Time Series") 


###############################################################################
# Spatial data

# spatial attributes automaticall identified
sosCoordinates(observations_nox)
# 1 line per feature, one feature per observation in the observation collection

# R spatial data structures: package sp
library("sp")
result_nox_sp <- as(observations_nox, "SpatialPointsDataFrame")
plot(result_nox_sp)

# 
bbox(coordinates(result_nox_sp))
#min       max
#x  30.96462  31.02729
#y -29.95696 -29.77789


###############################################################################
# Plot stations on a map
library(maps); library(mapdata); library(maptools); data(worldHiresMapEnv)
crs <- sosGetCRS(eo2h)[[1]]
worldHigh <- pruneMap(map(database = "worldHires",
				region = c("South Africa"),
#				ylim = c(-29.957, -29.776),
#				xlim = c(30.964, 31.029),
				plot = FALSE))
worldHigh_lines <- map2SpatialLines(worldHigh, proj4string = crs)

# plot:
plot(worldHigh_lines, col = "grey50",
		ylim = c(-29.97, -29.776),
		xlim = c(30.964, 31.029), axes = TRUE)

# add points to plot (stupid way, because multiple points for one location):
plot(result_nox_sp, add = TRUE, lwd = 2, pch = 21, col = "blue")
# unique(coordinates(result_nox_sp))

# add labels to plot:
labels_short <- lapply(sosCoordinates(observations_nox)[["feature"]],
		FUN = substring, first = 63)
text(x = sosCoordinates(observations_nox)[,1],
		y = sosCoordinates(observations_nox)[,2],
		pos = 4,
		labels = labels_short)


###############################################################################
# Plot stations on map with Google Maps background
library("raster"); library("dismo"); library("rgdal")

e <- extent(c(30.8, 31.2, -29.968, -29.77))
r <- raster(e, crs=proj4string(result_nox_sp))
values(r) <- runif(ncell(r))
g <- gmap(r, type = "terrain")
# ?Mercator
ptm <- spTransform(result_nox_sp,
		CRS("+proj=merc +a=6378137 +b=6378137 +lat_ts=0.0 +lon_0=0.0 +x_0=0.0 +y_0=0 +k=1.0 +units=m +nadgrids=@null +no_defs"))

# important when using StatET: Raster is not yet supported by that device
plot(g, interpolate = TRUE) #axes = TRUE)
points(ptm, col='red', pch="x", alpha=0.25, cex = 2)

ptmLabels <- Mercator(matrix(c(sosCoordinates(observations_nox)[,"lat"],
		sosCoordinates(observations_nox)[,"lon"]), ncol = 2))
text(x = ptmLabels[,"x"], y = ptmLabels[,"y"],
		pos = 4, labels = labels_short)


###############################################################################
# Climatology data for the state of Saxony, Germany
###############################################################################

# TODO analysis with dresden data!
sosOfferings(eo2h)["HUMIDITY"]
sosOfferings(eo2h)["AIR_PRESSURE"]
sosOfferings(eo2h)["TEMP"]

myTime <- sosCreateTime(sos = eo2h, time = "2003-01-01::2004-12-31")

# get data
observations_humidity <- getObservation(sos = eo2h,
		offering = "HUMIDITY",
		eventTime = myTime)
observations_air_pressure <- getObservation(sos = eo2h,
		offering = "AIR_PRESSURE",
		eventTime = myTime)
observations_temp <- getObservation(sos = eo2h,
		offering = "TEMP",
		eventTime = myTime)

# remove -999.0000 values, probably missing values
observations_temp_data <- sosResult(observations_temp)
hist(observations_temp_data[,"TEMP"])
observations_temp_data <- subset(observations_temp_data, TEMP > -900)
summary(observations_temp_data)
dim(sosResult(observations_temp)); dim(observations_temp_data)

observations_pressure_data <- sosResult(observations_air_pressure)
observations_pressure_data <- subset(observations_pressure_data, 
		AIR_PRESSURE > -900)
summary(observations_pressure_data)
hist(observations_pressure_data[,"AIR_PRESSURE"])
dim(sosResult(observations_air_pressure)); dim(observations_pressure_data)

observations_humidity_data <- sosResult(observations_humidity)
observations_humidity_data <- subset(observations_humidity_data, 
		HUMIDITY > -900)
summary(observations_humidity_data)
dim(sosResult(observations_humidity)); dim(observations_humidity_data)

# missing values
missing_values <- data.frame(
		c(dim(sosResult(observations_air_pressure))[[1]],
				dim(observations_pressure_data)[[1]]),
		c(dim(sosResult(observations_temp))[[1]],
				dim(observations_temp_data)[[1]]),
		c(dim(sosResult(observations_humidity))[[1]],
				dim(observations_humidity_data)[[1]]),
		row.names = c("all", "valid")
)
names(missing_values) <- c("air pressure", "temp", "humidity")
missing_values

# percentagea
missing_values["valid",] / missing_values["all",]
#		air pressure	temp		humidity
#valid	0.6767123		0.8599419	0.7695309

# barplot
barplot(as.matrix(missing_values), legend= rownames(missing_values),
		beside = TRUE)

# TODO continue with analysis

