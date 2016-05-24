# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")

################################################################################
# SOS @ CSIRO
# The South Esk test bed
cat("Go to the following website for details of the South Esk Hydrological Sensor Web - Tasmania, Australia: ",
		"http://www.csiro.au/sensorweb/au.csiro.OgcThinClient/OgcThinClient.html")

# Data subject to CSIRO's legal disclaimer: http://www.csiro.au/org/LegalNoticeAndDisclaimer.html

# See also:
# http://external.opengis.org/twiki_public/bin/view/ClimateChallenge2009/ServiceOfferingCSIRO
# www.csiro.au/sensorweb2/catalog/setup/ --> a catalog with the SOS urls!

csiro <- SOS("http://www.csiro.au/sensorweb/CSIRO_SOS/sos")

################################################################################
# TIME FORMAT ISSUES
################################################################################
#?Sys.timezone
Sys.timezone()
strptime("1995-05-25T15:30:00+11:00", format = "%Y-%m-%dT%H:%M:%OS")$isdst
# 1
strptime("1995-05-25T15:30:00+11:00", format = "%Y-%m-%dT%H:%M:%OS%z")$isdst
# -1
# isdst: Daylight Savings Time flag. Positive if in force, zero if not, negative if unknown.

#########################################
# This does not work on windows machines: the %z does not work for numerical
# outputs, so the create time string is never suitable and the timeFormat
# parameter cannot solve the problem.

# Try fixing with locale...
#Sys.getlocale(category = "LC_TIME")
#Sys.setlocale("LC_ALL", "English")
#Sys.getlocale()
#Sys.setenv(TZ="GMT") 
# All this does NOT remove the "Mitteleurop_ische Zeit" from the strftime output!!

################################################################################
# REPLACE TIME CONVERTER
csiroTimeConverter = function(x, sos) {
#	cat("Using adapted time parser for ", toString(x), "\n")
	.x <- paste (x, "00", sep = "")
	.value <- sosConvertTime(.x, sos = sos)
	return(.value)
}
# works:
#> csiroTimeConverter("2011-10-06T14:35:00+11", csiro)
#Using adapted time parser for  2011-10-06T14:35:00+11 
#[1] "2011-10-06 05:35:00 CEST"

# Ignore time zone when parsing, but use when creating output:
strftime(strptime("1995-05-25T15:30:00+11:00", format = sosDefaultTimeFormat), 
		format = paste(sosDefaultTimeFormat, "%z", sep  =""))
# Problem: Output is "1995-05-25T15:30:00Mitteleurop_ische Sommerzeit", not numerical!

##################
# REQUIRED OUTPUT (detected manually based on capabilities):
# 2011-10-06T10:28:10+11:00

################################################################################
# REPLACE only POSIXt encoder with a hack!
testtime <- strptime("2011-10-06T14:35:00+11", format = sosDefaultTimeFormat)
# returns POSIXlt - good.
strftime(testtime, format = sosDefaultTimeFormat)
format(testtime, format = sosDefaultTimeFormat)

setMethod(f = "encodeXML",
		signature = signature(obj = "POSIXt", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("CSIRO encoding... ")
			
			# time zone hack!
			.time <- obj + 11 * 60 * 60							# add 11 hours
			.formatted <- strftime(x = .time, format = sosTimeFormat(sos))
			.formatted <- paste(.formatted, "+11:00", sep = "")	# append 11:00
			
			if(verbose)
				cat("Formatted ", toString(obj), " to ", .formatted, "\n")
			
			return(.formatted)
		}
)

# CSIRO (orange, purple on map), replace testing one from above!
csiro <- SOS("http://www.csiro.au/sensorweb/CSIRO_SOS/sos",
		dataFieldConverters = SosDataFieldConvertingFunctions(
				"urn:ogc:data:time:iso8601" = csiroTimeConverter,
				"urn:ogc:def:phenomenon:OGC:rainfall" = sosConvertDouble),
		switchCoordinates = TRUE)

# Bureau of Meteorology (red and dark blue on map)
bom <- SOS("http://www.csiro.au/sensorweb/BOM_SOS/sos",
		dataFieldConverters = SosDataFieldConvertingFunctions(
				"urn:ogc:data:time:iso8601" = csiroTimeConverter,
				"urn:ogc:def:phenomenon:OGC:rainfall" = sosConvertDouble),
		switchCoordinates = TRUE, # verbose = TRUE
)
sosTimeFormat(bom)
sosOfferings(bom)[[1]]

# Tasmania Department of Primary Industries, Parks, Wildlife and Environment (DPIPWE, white on map)
dpiw <- SOS("http://www.csiro.au/sensorweb/DPIW_SOS/sos?Service=SOS&Request=GetCapabilities", switchCoordinates = TRUE)

# Hydro Tasmania Consulting - Remote Monitoring and Investigation Unit (yellow on map)
ht <- SOS("http://www.csiro.au/sensorweb/HT_SOS/sos")

# Forestry Tasmania - Fire Risk Management Branch (green on map)
forestry <- SOS("http://www.csiro.au/sensorweb/Forestry_SOS/sos")

# Tasmania Department of Primary Industries, Parks, Wildlife and Environment (DPIPWE) - Water Assessment Branch
dpipwe <- SOS("http://www.csiro.au/sensorweb/DPIPWE_SOS/sos")

# What about these?
#hutchins <- SOS("http://150.229.66.73/HutchinsSOS/sos")
#elliotwsn <- SOS("http://150.229.66.73/ElliotWSNSOS/sos")


################################################################################
# explore SOSs with plots

################################################################################
# plot two SOS with background map
crs <- sosGetCRS(csiro)

library(maps)
library(maptools)
world <- pruneMap(map(database = "world", region = "Australia:Tasmania",
				plot = FALSE))
world.lines <- map2SpatialLines(world, proj4string = crs[[2]])
plot(world.lines, col = "grey50")
plot(csiro, add = TRUE, lwd = 3, lty = 1)
plot(bom, add = TRUE,
		border.color.pal = sosDefaultColorPalette[[length(sosOfferings(csiro))]],
		lwd = 3, lty = 3)

# plot labels
labelCoords <- rbind(sosCoordinates(sosOfferings(csiro)),
		sosCoordinates(sosOfferings(bom)))
labels <- c(sosName(sosOfferings(csiro)), sosName(sosOfferings(bom)))
text(labels = labels, col = sosDefaultColorPalette,
		x = labelCoords[,1],
		y = labelCoords[,2])
title(main = paste(sosTitle(csiro), "and", sosTitle(bom)),
		sub = paste(sosAbstract(csiro), "\n", sosAbstract(bom)))

################################################################################
# plot one offering with high resolution background map and cities, including 
# map axes and scale
off <- sosOfferings(csiro)[[2]]
library(mapdata)
data(worldHiresMapEnv)

# detect region automatically
region <- map.where(database = "worldHires", sosCoordinates(off))
worldHigh <- pruneMap(map(database = "worldHires", region = region,
				plot = FALSE))
worldHigh_Lines <- map2SpatialLines(worldHigh, proj4string = crs[[2]])

plot(worldHigh_Lines, col = "grey50")
data(world.cities)
map.cities(label = TRUE, pch = 19, col = "black")

plot(off, add = TRUE, border = "red", lwd = 3)
title(main = paste("Offering with ID '", sosId(off), "'", sep = ""),
		sub = paste("Features:", toString(sosFeaturesOfInterest(off))))

map.axes()
map.scale(metric = TRUE, ratio = FALSE)


################################################################################
# Data request and consolidation based on three SOS: bom, csiro and dpiw

# phenomenon rainfall or rainfalltoday is available at all stations
rainfall <- "urn:ogc:def:phenomenon:OGC:rainfall"

lastDay <- sosCreateTimePeriod(sos = bom, begin = (Sys.time() - 3600 * 24),
		end = Sys.time())
sosTimeFormat(bom); encodeXML(lastDay, bom)
str(lastDay)

#####
# bom
#
rainfall.off.bom <- sosOfferings(bom)[["BOM Offering"]]
phenomenon.bom <- as.list(unlist(sosObservedProperties(bom)))
phenomenon.bom <- phenomenon.bom[grep(pattern = "rain", phenomenon.bom)]

rainfall.obs.bom <- getObservation(sos = bom, offering = rainfall.off.bom, 
		observedProperty = phenomenon.bom, verbose = TRUE,
		eventTime = sosCreateEventTimeList(lastDay))
rainfall.result.bom <- sosResult(rainfall.obs.bom, coordinate = TRUE)
summary(rainfall.result.bom)

########################
# csiro, ht and datacell
#
rainfall.off.csiro <- sosOfferings(csiro)[["Rain Gauges"]]
rainfall.off.ht <- sosOfferings(csiro)[["HT Weather Stations"]]
rainfall.off.datacell <- sosOfferings(csiro)[["Datacall Weather Stations"]]

procedures.csiro <- sosProcedures(rainfall.off.csiro)
phenomenon.csiro <- sosObservedProperties(rainfall.off.csiro)
phenomenon.csiro <- phenomenon.bom[grep(pattern = "rain", phenomenon.csiro)]
print(paste("FOIs: ", toString(sosFeaturesOfInterest(rainfall.off.csiro))))

rainfall.obs.csiro <- getObservation(sos = csiro, # verbose = TRUE,
		offering = rainfall.off.csiro,
		procedure = procedures.csiro, 
		observedProperty = phenomenon.csiro,
		eventTime = sosCreateEventTimeList(lastDay))
# at this point realized problem with parsing time (just NAs), need to fix it
# with csiroTimeParser (see above)

phenomenon.ht <- sosObservedProperties(rainfall.off.ht)
phenomenon.ht <- phenomenon.ht[grep(pattern = "rain", phenomenon.ht)]
rainfall.obs.ht <- getObservation(sos = csiro, # verbose = TRUE
		offering = rainfall.off.ht, 
		observedProperty = phenomenon.ht,
		eventTime = sosCreateEventTimeList(lastDay))

phenomenon.dc <- sosObservedProperties(rainfall.off.datacell)
phenomenon.dc <- phenomenon.dc[grep(pattern = "rainfalltoday", phenomenon.dc)]
rainfall.obs.dc <- getObservation(sos = csiro, # verbose = TRUE
		offering = rainfall.off.datacell, 
		observedProperty = phenomenon.dc,
		eventTime = sosCreateEventTimeList(lastDay))

# get data values
rainfall.result.csiro <- sosResult(rainfall.obs.csiro, coordinate = TRUE)
summary(rainfall.result.csiro)

rainfall.result.ht <- sosResult(rainfall.obs.ht, coordinate = TRUE)
summary(rainfall.result.ht)

#crs.ht <- sosGetCRS(rainfall.obs.ht[[1]])
rainfall.result.dc <- sosResult(rainfall.obs.dc, coordinate = TRUE)
summary(rainfall.result.dc)

################################################################################
# Data consolidation:
# save all data in analyzable data structure for one point in time.

# Times do not match exactly:
#time.csiro <- rainfall.result.csiro[,"Time"]
#time.ht <- rainfall.result.ht[,"Time"]
#time.dc <- rainfall.result.dc[,"Time"]
#time.bom <- rainfall.result.bom[,"Time"]

#
# subset data for one point in time
#
rainfall.data.time.bom <- unique(rainfall.result.bom[,"Time"])[[1]]
rainfall.data.bom <- subset(rainfall.result.bom,
		Time == rainfall.data.time.bom,
		c("lat", "lon", "feature", "urn:ogc:def:phenomenon:OGC:rainfall"))
rainfall.data.bom <- cbind(rainfall.data.bom, offering = c("bom"))

rainfall.data.time.csiro <- unique(rainfall.result.csiro[,"Time"])[[1]]
rainfall.data.csiro <- subset(rainfall.result.csiro,
		Time == rainfall.data.time.csiro,
		c("lat", "lon", "feature", "urn:ogc:def:phenomenon:OGC:rainfall"))
rainfall.data.csiro <- cbind(rainfall.data.csiro, offering = c("csiro"))

# some error here, remove for now
#rainfall.data.time.ht <- unique(rainfall.result.ht[,"Time"])[[1]]
#rainfall.data.ht <- subset(rainfall.result.ht,
#		Time == rainfall.data.time.ht,
#		c("lat", "lon", "feature", "urn:ogc:def:phenomenon:OGC:rainfall"))
#rainfall.data.ht <- cbind(rainfall.data.ht, offering = c("ht"))

rainfall.data.time.dc <- unique(rainfall.result.dc[["Time"]])[[1]]
rainfall.data.dc <- subset(rainfall.result.dc,
		Time %in% rainfall.data.time.dc,
		c("lat", "lon", "feature", "urn:ogc:def:phenomenon:OGC:rainfalltoday"))
names(rainfall.data.dc) <- c(names(rainfall.data.dc)[c(1,2,3)],
		"urn:ogc:def:phenomenon:OGC:rainfall")
rainfall.data.dc <- cbind(rainfall.data.dc, offering = c("dc"))

#
# Bind data frames together
#
rainfall.data <- rbind(rainfall.data.bom, 
		rainfall.data.csiro,
#		rainfall.data.ht,
		rainfall.data.dc)
names(rainfall.data) <- list("lat", "lon", "feature", "rainfall", "offering")
summary(rainfall.data)
# see offering!

# Create SpatialPointsDataFrame, use just the CRS of one observation
crs <- sosGetCRS(rainfall.obs.csiro)[[2]]
rainfall.spdf <- SpatialPointsDataFrame(
		coords = rainfall.data[,c("lon", "lat")],
		data = rainfall.data[,c("feature", "rainfall", "offering")],
		proj4string = crs)
str(rainfall.spdf)
bbox(rainfall.spdf)
summary(rainfall.spdf)

cat("Please be aware that the times that are combined in the following are NOT equal!\n")
cat(paste("BOM:\t", toString(rainfall.data.time.bom), "\n"),
		paste("CSIRO:\t", toString(rainfall.data.time.csiro), "\n"),
#		paste("HT:\t", toString(rainfall.data.time.ht), "\n"),
		paste("DC:\t", toString(rainfall.data.time.dc), "\n"))


################################################################################
# plot stations with background data
library("maps")
library("sp")
library("rgdal") # for spTransform
library("maptools") # for pruneMap
library("mapdata")

world.p <- pruneMap(
		map(database = "worldHires", region = "Australia:Tasmania",
				plot = FALSE))
world.sp <- map2SpatialLines(world.p, proj4string = crs)
plot(x = world.sp, col = "grey", main = "Rainfall Tasmania")
plot(rainfall.spdf, add = TRUE)

#text(x = coordinates(rainfall.spdf)[,"lat"],
#		y = coordinates(rainfall.spdf)[,"lon"],
#		labels = rainfall.spdf@data[, "feature"], adj = c(0, 1), cex = 0.75)
text(x = coordinates(rainfall.spdf)[,"lon"],
		y = coordinates(rainfall.spdf)[,"lat"], col = "blue",
		labels = rainfall.spdf@data[,"offering"], adj = c(1, 0), cex = 0.75)
title(main = paste("Observation locations of", sosTitle(csiro)),
		sub = paste(sosAbstract(csiro)))

# inspect and plot data
summary(rainfall.spdf[,"rainfall"])
hist(rainfall.spdf@data[,"rainfall"], main = "Rainfall")
bubble(rainfall.spdf, zcol = 2, maxsize = 2, col = c("#ff5588", "#ff5588"),
		main = "Rainfall in Tasmania", do.sqrt = TRUE)

################################################################################
# transform to UTM for kriging and background map
utm55 = CRS("+proj=utm +zone=55 +datum=WGS84")
bgmap = map2SpatialLines(map("worldHires", region = "Australia:Tasmania",
				plot=F))
proj4string(bgmap) <- "+proj=longlat +datum=WGS84"
bgmap.utm <- spTransform(bgmap, utm55)
rainfall.spdf.utm = spTransform(x = rainfall.spdf, CRSobj = utm55)
# error with spTransform, it makes all coordinates the same!
plot(bgmap.utm, col = "grey")
plot(rainfall.spdf.utm, add = TRUE)

################################################################################
# Intamap
#library("intamap")
#one.intamap <- one.spdf
#names(one.intamap) <- list("feature", "value")
#
#obj <- createIntamapObject(observations = one.intamap,
#		targetCRS = crs.csiro)
#checkSetup(obj)
# less than 20 observations, not able to perfrom interpolation

################################################################################
# automap
library("automap")
# Ordinary kriging, no new_data object
kriging_result = autoKrige(
		log(rainfall.spdf[["rainfall"]])~1,
		input_data = rainfall.spdf,
		new_data = SpatialPixels(
				SpatialPoints(makegrid(rainfall.spdf, n = 300))))
plot(kriging_result)

################################################################################
# Kriging
library("gstat")
rainfall.grid.utm = SpatialPixels(
		SpatialPoints(makegrid(rainfall.spdf.utm, n = 300), 
				proj4string = utm55))
m <- vgm(.59, "Sph", 874, .04)
# ordinary kriging:
x <- krige(log(rainfall.spdf[["rainfall"]])~1,
		rainfall.spdf.utm, rainfall.grid.utm, model = m)
spplot(x["var1.pred"], main = "ordinary kriging predictions")
spplot(one.utm, add = TRUE)

spplot(x["var1.var"],  main = "ordinary kriging variance")

################################################################################
# spacetime
# TODO continue analysis with spacetime package or make a forecast 
# See http://robjhyndman.com/researchtips/forecast3/

###################################
# Demo finished, try another one! #
###################################
