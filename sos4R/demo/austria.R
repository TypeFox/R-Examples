# Copyright (C) 2011 by 52 North Initiative for Geospatial Open Source Software GmbH, Contact: info@52north.org
# This program is free software; you can redistribute and/or modify it under the terms of the GNU General Public License version 2 as published by the Free Software Foundation. This program is distributed WITHOUT ANY WARRANTY; even without the implied WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details. You should have received a copy of the GNU General Public License along with this program (see gpl-2.0.txt). If not, write to the Free Software Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or visit the Free Software Foundation web page, http://www.fsf.org.
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r
library("sos4R")

################################################################################
# http://ispacevm10.researchstudio.at/sostester/
# Contact: Michael Lippautz: michael.lippautz@researchstudio.at
# Use the two development services on ispacevm10.

# PROBLEMS using GET:
# - Temporal filtering has no effect...

################################################################################
# Nationalpark Berchtesgaden
npbg_converter <- SosDataFieldConvertingFunctions(
		"urn:ogc:def:property:OGC:Time:iso8601" = sosConvertTime,
		"urn:ogc:def:property:OGC:Reflection" = sosConvertDouble,
		"urn:ogc:def:property:OGC:Insolation" = sosConvertDouble)
npbg <- SOS("http://ispacevm10.researchstudio.at/geoservices/npbg",
		method = "GET",
		#verboseOutput = TRUE,
		dataFieldConverters = npbg_converter,
		sections = NA)
npbg
summary(npbg)

#################
# Plot whole sos:
library(maps); library(mapdata); library(maptools)
data(worldHiresMapEnv)
crs <- unique(sosGetCRS(npbg))[[1]]
worldHigh <- pruneMap(map(database = "worldHires", region = "Austria",
				plot = FALSE))
worldHigh_Lines <- map2SpatialLines(worldHigh, proj4string = crs)

plot(worldHigh_Lines, col = "grey50")
plot(npbg, add = TRUE, lwd = 2)
map.axes()
map.scale()
offNames <- sapply(names(sosOfferings(npbg)), FUN = strsplit, split = ":")
title(main = paste("Offerings by '", sosTitle(npbg), "'", sep = ""),
		sub = toString(sapply(offNames, "[[", 3)))

##########
# zoom in:
poly1 <- as(sosOfferings(npbg)[[1]], "Spatial")
bbox(poly1)

bounds <- sosBoundedBy(sosOfferings(npbg), bbox = TRUE)
bounds_xlim <- c(min(sapply(bounds, "[[", "coords.lon", "min")) - 0.2,
		max(sapply(bounds, "[[", "coords.lon", "max")) + 0.2)
bounds_ylim <- c(min(sapply(bounds, "[[", "coords.lat", "min") - 0.2),
		max(sapply(bounds, "[[", "coords.lat", "max")) + 0.2)
# increase a little bit (works only because all > 0

# run same code again:
plot(worldHigh_Lines, col = "grey50", xlim = bounds_xlim, ylim = bounds_ylim)
plot(npbg, add = TRUE, lwd = 2)
map.axes()
map.scale()
offNames <- sapply(names(sosOfferings(npbg)), FUN = strsplit, split = ":")
title(main = paste("Offerings by '", sosTitle(npbg), "'", sep = ""),
		sub = toString(sapply(offNames, "[[", 3)))

# plot stations NOT POSSIBLE because SOS expects parameter "SensorId", not
# "procedure":
# http://ispacevm10.researchstudio.at/geoservices/npbg?service=SOS&request=DescribeSensor&version=1.0.0&SensorId=org%3Anpbg%3ABlaueis&outputFormat=text/xml;subtype%3D%22sensorML/1.0.1%22
#describeSensor(npbg, procedure = "org:npbg:Blaueis", verbose = TRUE)
#npbg_procedures <- unique(unlist(sosProcedures(npbg)))
#procs_descr <- lapply(X = npbg_procedures, FUN = describeSensor, # verbose = TRUE,
#		sos = npbg)
#procs_descr[[1]]
#for (x in procs_descr) {
#	plot(x, add = TRUE, pch = 19)
#}
#text(sosCoordinates(procs_descr)[c("x", "y")], labels = sosId(procs_descr),
#		pos = 4)


#########################
# superordinate offering:
np.off <- sosOfferings(npbg)[["org:npbg:Nationalpark"]]
#np.off
summary(np.off)

npbg_obsProp <- sosObservedProperties(np.off)
npbg_obsProp
npbg_proc <- sosProcedures(np.off)
npbg_proc

###########
# Get data:
lastDay <- sosCreateTimePeriod(sos = npbg, begin = (Sys.time() - 3600 * 24),
		end = Sys.time())

obs_proc1 <- getObservation(sos = npbg, offering = np.off, inspect = TRUE,
	procedure = npbg_proc[[1]],
	eventTime = sosCreateEventTimeList(lastDay)
	)

# cannot get the coordinates automatically via sosResult(), because they are
# given as attribute, not as a featureOfInterest
#result.proc1 <- sosResult(obs.proc1, coordinates = TRUE)
#spdf.proc1 <- as(obs.proc1, "Spatial")

# Coordinates inline:
result_proc1 <- sosResult(obs_proc1)
summary(result_proc1)

coords.proc1 <- unique(result_proc1[c("Latitude", "Longitude")])
coords.proc1


##################
# plot all values:
names(result_proc1)
plot(result_proc1[8:9],
		main = paste(npbg_proc[[1]], "at", toString(coords.proc1), "from, to",
				toString(range(result_proc1[["Time"]]))))

######################
# xyplot, dotplot ...:
xyplot(Insolation ~ RelativeHumidity, data = result_proc1, 
		main = paste(npbg_proc[[1]], "(", 
				toString(coords.proc1), ")"),
		sub = paste("Time range: ", toString(range(result_proc1[["Time"]]))))

xyplot(Insolation ~ RelativeHumidity | AirTemperature,
		data = result_proc1[1:100,])

dotplot(Insolation ~ RelativeHumidity, data = result_proc1)


####################################
# plot values against time with xts:
result_proc1_times <- unique(result_proc1[["Time"]])

library(xts)
proc1_xts_refl <- xts(x = result_proc1[["Reflection"]],
		order.by = result_proc1[["Time"]], unique = TRUE)
proc1_xts_inso <- xts(x = result_proc1[["Insolation"]],
		order.by = result_proc1[["Time"]], unique = TRUE)
plot(proc1_xts_refl,
		main = paste("Reflection at", npbg_proc[[1]], "(", 
				toString(coords.proc1), ")"),
		type = "bars")

############################
# plot time series with zoo:
library(zoo)

proc1_zoo <- zoo(x = as.matrix(result_proc1[5:9]),
		order.by = result_proc1[["Time"]])
str(proc1_zoo)
plot(proc1_zoo, main = paste("Time Series at", npbg_proc[[1]], "(", 
				toString(coords.proc1), ")"),
		plot.type = "multiple")


####################
# request more data:
obs_proc123 <- getObservation(sos = npbg, offering = np.off, # inspect = TRUE,
		procedure = npbg_proc[1:3])
#Finished getObservation to http://ispacevm10.researchstudio.at/geoservices/npbg 
#--> received 3 observation(s) having 31345 result values [ 8942, 11112, 11291 ].
obs_proc123
str(obs_proc123[[1]])
sosProcedures(obs_proc123)

################################################################################
# Land Oberoesterreich
ooe <- SOS("http://ispacevm10.researchstudio.at/geoservices/ooe", sections = NA,
		verboseOutput = TRUE)
ooe
summary(ooe)

# TODO KML export, or make visualization example following howto:
# http://spatial-analyst.net/wiki/index.php?title=Export_maps_to_GE


###################################
# Demo finished, try another one! #
###################################
