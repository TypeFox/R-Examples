# Copyright (C) 2010 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: 	Daniel Nuest (daniel.nuest@uni-muenster.de)
#			Edzer Pebesma (edzer.pebesma@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")

################################################################################
# AQE SOS
# Data source: http://www.eea.europa.eu/themes/air/airbase

# Set the converters for observed properties:
aqe.converters <- SosDataFieldConvertingFunctions(
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[PM10]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[NO2]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[O3]" = sosConvertDouble,
		"http://www.opengis.net/def/property/OGC/0/SamplingTime" = sosConvertTime,
		"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString)

# Create the SOS connection:
aqe <- SOS(url = "http://giv-uw.uni-muenster.de:8080/AQE/sos",
		dataFieldConverters = aqe.converters)
summary(aqe)

# Get the available offerings:
aqe_offerings <- sosOfferings(aqe)
names(aqe_offerings)

###########
# Plot SOS:
library(maps); library(mapdata); library(maptools)
if(!require(rgdal, quietly = TRUE))
	print("rgdal not present: CRS values will not be converted correctly")
data(worldHiresMapEnv)
crs <- sosGetCRS(aqe)[[1]]
region <- map.where(database = "worldHires",
		sosCoordinates(aqe_offerings)) # find region
worldHigh <- pruneMap(map(database = "worldHires", region = region,
				plot = FALSE))
worldHigh_Lines <- map2SpatialLines(worldHigh, proj4string = crs)

plot(worldHigh_Lines, col = "grey50")
plot(aqe, add = TRUE, lwd = 3)
title(main = paste("Offerings Germany by '", sosTitle(aqe), "'", sep = ""),
		sub = toString(names(aqe_offerings)))

################################################################################
# NO2
#
# Extract one offering of interest and explore:
aqe_off_no2 <- aqe_offerings[["NO2"]]
aqe_off_no2
summary(aqe_off_no2)

# Get observations, december 2003 is arbitrary choice!
dec2003.12Hrs = sosCreateTime(sos = aqe,
		time = "2003/12/01 08:00::2003/12/01 20:00")
dec2003.24Hrs = sosCreateTime(sos = aqe,
		time = "2003/12/01 08:00::2003/12/02 08:00")
dec2003 = sosCreateTime(sos = aqe, time = "2003/12/01::2003/12/31")

# Request data (request and response can be check by setting the inspect=TRUE):
obs_no2_12Hrs <- getObservation(sos = aqe, # inspect = TRUE,
		offering = aqe_off_no2,
		#procedure = sosProcedures(aqe.off.no2)[1:20],
		saveOriginal = TRUE, # saves file in getwd()
		eventTime = dec2003.12Hrs)
# 38 secs

################################################
# Reloading the file later for further analysis:
# Use the file name printed out by the getObservation(...) call or check your
# working directory:
list.files(path = getwd(), pattern = sosName(aqe_off_no2))
#parseFile("NO2_2011-03-03_11:31:23.xml")

##############################################
# Explore the returned observation collection:
obs_no2_12Hrs
# There is one observatio for every FOI / procedure combination :
names(obs_no2_12Hrs)[1:3]

########################
# Subset the collection:
# (features, observed properties and procedures):
# sosFeatureIds(obs.no2.12Hrs)[c(1,100)]
obs_no2_12Hrs[sosFeatureIds(obs_no2_12Hrs)[c(1,100)]]
# sosObservedProperties(obs.no2.12Hrs)[2]
obs_no2_12Hrs[sosObservedProperties(obs_no2_12Hrs)[2]]
# sosProcedures(obs.no2.12Hrs)[200:201]
obs_no2_12Hrs[sosProcedures(obs_no2_12Hrs)[200:201]]

# More requests for more data for testing of response time:
#obs.no2.24Hrs <- getObservation(sos = aqe, # inspect = TRUE,
#		offering = aqe.off.no2,
#		eventTime = dec2003.24Hrs)
## 41 secs
#obs.no2.dec <- getObservation(sos = aqe, # inspect = TRUE,
#		offering = aqe.off.no2,
#		eventTime = dec2003)
## 3:25 mins
#obs.no2.dec
#result.no2.dec <- sosResult(obs.no2.dec)[1:10, ]


################################################################################
# Get the result data for all observations, with coordinates:
result_no2_12Hrs <- sosResult(obs_no2_12Hrs, coordinates = TRUE)
# Coordinates only:
#sosCoordinates(obs.no2.12Hrs[1:10])
# One observation only
# sosResult(obs.no2.12Hrs[[42]])

summary(result_no2_12Hrs)
NO2 <- colnames(result_no2_12Hrs)[[3]]

#################################################
# Subset and sort the data with subset or sort_df
subset(result_no2_12Hrs, feature=="foi_DEBY109")
library("reshape")
# The ten highest values:
tail(sort_df(result_no2_12Hrs, NO2), 10)

########################
# Histogram of NO2 data:
hist(result_no2_12Hrs[,3], main = "NO2")
# Test plot:
plot(result_no2_12Hrs[["SamplingTime"]], result_no2_12Hrs[[NO2]])


################################################################################
# Get the result data and create sp object:
obs.no2.crs <- sosGetCRS(obs_no2_12Hrs)
no2_spdf <- SpatialPointsDataFrame(
		coords = result_no2_12Hrs[,c("lon", "lat")],
		data = result_no2_12Hrs[,c("SamplingTime", "feature", NO2)],
		proj4string = obs.no2.crs)
bbox(no2_spdf)
#obs.no2.bbox <- sosBoundedBy(obs.no2.12Hrs, bbox = TRUE) # equal
summary(no2_spdf)

#########################################
# Shortcut to get SpatialPointsDataFrame:
#as(obs.no2.12Hrs[[1]], "SpatialPointsDataFrame")
no2_spdf_shortcut <- as(obs_no2_12Hrs, "SpatialPointsDataFrame")
summary(no2_spdf_shortcut)

####################################
# Plot stations with background map:
library("mapdata")
germany_p <- pruneMap(map(database = "worldHires", region = "Germany",
				plot = FALSE))
germany_sp <- map2SpatialLines(germany_p, proj4string = obs.no2.crs)
proj4string(germany_sp) <- obs.no2.crs
plot(x = germany_sp, col = "grey")
plot(no2_spdf, pch = 20, col = "blue", add = TRUE)
title("NO2 Germany")

##############
# Bubble plot:
bubble(no2_spdf, zcol = 3, maxsize = 2, col = c("#1155ff"),
		main = "NO2 in Germany", do.sqrt = TRUE)

################################################################################
# Transform to UTM for kriging and background map:
utm32 = CRS("+proj=utm +zone=32 +datum=WGS84")
germany_utm <- spTransform(germany_sp, utm32)
no2_spdf_utm = spTransform(no2_spdf, utm32)
plot(germany_utm, col = "grey")
points(no2_spdf_utm, cex=.5, pch=3)
title(main = "NO2 Sensor Stations, Germany", sub = "UTM projection")
 
################################################################################
# Spatial interpolation, averaging over (essentially ignoring) time:
library(cshapes) # get a polygon of Germany, rather than a set of lines:
cs = cshp()
g = spTransform(cs[cs$CNTRY_NAME=="Germany",], utm32)
# create a grid of points within Germany:
grdpoints = SpatialPoints(makegrid(germany_utm),utm32)
grd = SpatialPixels(grdpoints)[g]
names(no2_spdf_utm)[3] = "NO2"
library(gstat)
no2.id = idw(NO2~1, no2_spdf_utm, grd)
lt=list(list("sp.polygons", g),
	list("sp.points", no2_spdf_utm, col=grey(.7), cex=.5))
spplot(no2.id[1], sp.layout = lt, col.regions = bpy.colors(),
		main = paste("IDW Interpolation of NO2,",
				min(no2_spdf_utm$SamplingTime), "to",
				max(no2_spdf_utm$SamplingTime)))

## now get the values of 2003-12-01 12:00:00 CET, or the third time stamp
# by creating a spatio-temporal structure, and indexing the time axis:
t12h = unique(no2_spdf_utm$SamplingTime)[3]
library(spacetime)
no2_stidf = STIDF(as(no2_spdf_utm, "SpatialPoints"), 
	no2_spdf_utm$SamplingTime, data.frame(NO2 = no2_spdf_utm$NO2))
no2_stfdf = as(no2_stidf, "STFDF")
no2_12h = no2_stfdf[,3,"NO2"]
no2_12h2 = no2_stfdf[,t12h,"NO2"]
all.equal(no2_12h, no2_12h2)
# inverse distance interpolation of the 12:00h NO2 values:
no2_12h_id = idw(NO2~1, no2_12h[!is.na(no2_12h$NO2),], grd)
spplot(no2_12h_id[1], sp.layout = lt, col.regions = bpy.colors(),
		main = paste("IDW Interpolation of NO2,", t12h))

################################################################################
# Plot with whole year 2004 for one station:
# See http://www.eea.europa.eu/themes/air/airbase/map-stations.
denw095 <- "urn:ogc:object:feature:Sensor:EEA:airbase:4.0:DENW095"
denw095_descr <- describeSensor(aqe, denw095)
denw095_descr
#denw095.descr@xml

# Get the identifier of the station:
denw095_id <- xmlValue(getNodeSet(doc = denw095_descr@xml,
		path = "//sml:Term[@definition='urn:ogc:def:identifier:OGC:1.0:longName']/sml:value/text()",
		namespaces = sos4R:::.sosNamespaceDefinitionsSML)[[1]])

# Request observations:
obs_denw095_2004 <- getObservation(sos = aqe, # inspect = TRUE,
		offering = aqe_off_no2,
		procedure = denw095,
		eventTime = sosCreateEventTimeList(sosCreateTimePeriod(sos = aqe,
						begin = as.POSIXct("2004/01/01"),
						end = as.POSIXct("2004/12/31")))
)

# Print statistical information and plot time series:
data_denw095_2004 <- sosResult(obs_denw095_2004)
summary(data_denw095_2004)

denw095.NO2.attributes <- attributes(data_denw095_2004[[NO2]])
plot(data_denw095_2004[["SamplingTime"]], data_denw095_2004[[NO2]], type = "l",
		main = paste("NO2 in", denw095_id, "2004"), sub = denw095,
		xlab = "Time",
		ylab = paste("NO2 (",
				denw095.NO2.attributes[["unit of measurement"]],
				")", sep = ""))
data.denw095.2004.locRegr = loess(data_denw095_2004[[NO2]]~as.numeric(data_denw095_2004[["SamplingTime"]]),
		data_denw095_2004, enp.target = 30)
p = predict(data.denw095.2004.locRegr)
lines(p ~ data_denw095_2004[["SamplingTime"]], col = 'blue',lwd = 4)

###################################
# Demo finished, try another one! #
###################################
