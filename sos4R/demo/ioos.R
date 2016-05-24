# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")


################################################################################
# IOOS/NDBC
#
# Service Description: http://sdf.ndbc.noaa.gov/sos/
# Test Client: http://sdf.ndbc.noaa.gov/sos/test.shtml
#
# IOOS Map: http://www.ioos.gov/catalog/
# NDBC Map: http://www.ndbc.noaa.gov/
# Good for selecting subgroups/programmes of sensors
#
ioos <- SOS(url = "http://sdf.ndbc.noaa.gov/sos/server.php",
		timeFormat = "%Y-%m-%dT%H:%M:%SZ")
ioos.off <- sosOfferings(ioos)
names(ioos.off)
sosResponseFormats(ioos)

#############
# first test:
# !!! use sosName() of offering
getObservation(ioos, offering = sosName(ioos.off[[100]]), # verbose = TRUE
		responseFormat = "text/csv",
		observedProperty = sosObservedProperties(ioos.off[[100]])[1])

off.coords <- sosCoordinates(ioos.off)
off.names <- sosName(ioos.off)

# Plot offerings:
crs <- sosGetCRS(ioos.off[[1]])
library(maps); library(mapdata); library(maptools)
world <- pruneMap(map(database = "world", plot = FALSE))
world.lines <- map2SpatialLines(world, proj4string = crs)

plot(world.lines, col = "grey50")
plot(ioos, lwd = 3, add = TRUE)
title(main = sosTitle(ioos))
#text(x = off.coords[,1], y = off.coords[,1], col = "black",
#		labels = off.names, adj = c(1, 0), cex = 0.75)

ioos.procedures <- unique(unlist(sosProcedures(ioos)))
length(ioos.procedures); length(ioos.off)
# seems like one offering per sensor

################################################################################
# get data:

# offerings in pacific?
# TODO all also the ones starting with 51, 41, 32, ...

# TODO get all TAO/TRITON Array sensors: http://www.pmel.noaa.gov/tao/index.shtml
# Terminology: http://www.pmel.noaa.gov/tao/jsdisplay/help/help_terminology_f.html
# Data display: http://www.pmel.noaa.gov/tao/jsdisplay/
# Create similar plot: http://www.pmel.noaa.gov/cgi-tao/cover.cgi?P1=uwnd&P2=20110304-March-6-2011&P3=month&P4=off&script=jsdisplay/scripts/lat-lon-5day-jsd.csh
# About the buoys: http://www.pmel.noaa.gov/tao/proj_over/pubs/mil96paper.html

# OR JUST GET LATEST VALUE FOR ALL INSTEAD!!
offerings.wmo52 <- ioos.off[grep(pattern = "wmo:52", x = off.names)]
obsProps.wmo52 <- unique(unlist(sosObservedProperties(offerings.wmo52)))

last48hrs = sosCreateEventTimeList(sosCreateTimePeriod(ioos,
				begin = as.POSIXct(Sys.time() - 3600 * 48),
				end = as.POSIXct(Sys.time())))
phenomenon <- list("http://mmisw.org/ont/cf/parameter/sea_water_temperature")

# use lapply to call getObservation for every offering with "...:wmo:52..."
# possible alternative: use "all" offering:
#sosName(ioos.off[[1]])
obs.wmo52 <- lapply(X = sosName(offerings.wmo52), FUN = getObservation, #verbose = TRUE, 
		sos = ioos, observedProperty = phenomenon, 
		eventTime = last48hrs, 
		responseFormat = "text/csv")
length(obs.wmo52)
obs.wmo52[[1]]
#attributes(obs.wmo52[[1]])
names(obs.wmo52[[1]])

obs.wmo52.all <- sosResult(obs.wmo52)
summary(obs.wmo52.all)
str(obs.wmo52.all)
# columns are all factors, convert!
obs.wmo52.all[["sea_water_temperature (C)"]] <- as.numeric(obs.wmo52.all[["sea_water_temperature (C)"]])
obs.wmo52.all[["latitude (degree)"]] <- as.numeric(obs.wmo52.all[["latitude (degree)"]])
obs.wmo52.all[["longitude (degree)"]] <- as.numeric(obs.wmo52.all[["longitude (degree)"]])
obs.wmo52.all[["depth (m)"]] <- as.numeric(obs.wmo52.all[["depth (m)"]])
obs.wmo52.all[["date_time"]] <- as.POSIXct(obs.wmo52.all[["date_time"]])
summary(obs.wmo52.all)
dim(obs.wmo52.all)

hist(obs.wmo52.all[["sea_water_temperature (C)"]])
colnames(obs.wmo52.all)

# coordinates seem wrong and there are NA values in the 
# data.frame, must be removed
obs.wmo52.all <- obs.wmo52.all[complete.cases(obs.wmo52.all),]

spdf <- SpatialPointsDataFrame(
		coords = obs.wmo52.all[c(4,3)],
		data = obs.wmo52.all[-c(4,3)],
		proj4string = crs)
summary(spdf)

plot(world.lines, col = "grey50")
plot(spdf, add = TRUE)


################################################################################
# most recent observation:
obs.csv <- getObservation(ioos, offering = sosName(ioos.off[[100]]),
		responseFormat = "text/csv",
		observedProperty = sosObservedProperties(ioos.off[[100]])[2])
obs.csv
summary(obs.csv)


################################################################################
# describe sensor:
# requires SensorML 1.0.0
ioos.get <- SOS(url = "http://sdf.ndbc.noaa.gov/sos/server.php",
		method = SosSupportedConnectionMethods()[["GET"]],
		timeFormat = "%Y-%m-%dT%H:%M:%SZ")
describeSensorOp <- sosOperation(ioos.get, sosDescribeSensorName)
describeSensor.outputFormat <- describeSensorOp@parameters[["outputFormat"]][[1]]
ioos.procedures <- unique(unlist(sosProcedures(ioos.get)))
ioos.sensor.1.1 <- describeSensor(sos = ioos.get, procedure = ioos.procedures[[1]][[1]],
		outputFormat = describeSensor.outputFormat, verbose = TRUE)
ioos.sensor.1.1
ioos.sensor.1.1@xml

##############
# time series:
#begin <- sosTime(ioos.post.off[[1]], convert = TRUE)[[1]]
end <- as.POSIXct(Sys.time())
begin <- end - 3600 * 24 * 30

.offering <- sosOfferings(ioos, name = sosName(ioos.off[[1]]))
.offeringId <- sosId(.offering)
obs.001 <- getObservation(sos = ioos,
		offering = sosName(.offering), # "urn:ioos:network:noaa.nws.ndbc:all"
		procedure = sosProcedures(ioos.post.off[[1]])[690:700],
		observedProperty = sosObservedProperties(ioos.post.off[[1]])[6:7],
		responseFormat = "text/csv",
		eventTime = sosCreateEventTimeList(
				sosCreateTimePeriod(sos = ioos.post, begin = begin, end = end)),
		inspect = TRUE, verbose = TRUE)


################################################################################
# KML
kml <- getObservation(ioos, offering = "urn:ioos:network:noaa.nws.ndbc:all",
#		verbose = TRUE,
#		saveOriginal = TRUE,
		responseFormat = "application/vnd.google-earth.kml+xml",
		observedProperty = list(
				"http://mmisw.org/ont/cf/parameter/air_temperature"))
kml

# TODO do sth. with the KML, e.g. export using examples from Spatial-Analyst?
# TODO use plotKML to create an output?


################################################################################
# GET
ioos.get <- SOS(url = "http://sdf.ndbc.noaa.gov/sos/server.php",
		method = SosSupportedConnectionMethods()[["GET"]],
		timeFormat = "%Y-%m-%dT%H:%M:%SZ")
#		parsers = SosParsingFunctions("GetObservation" = parseNoParsing)


###################################
# Demo finished, try another one! #
###################################
