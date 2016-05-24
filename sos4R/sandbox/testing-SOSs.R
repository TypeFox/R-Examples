################################################################################
# Copyright (C) 2010 by 52 North                                               #
# Initiative for Geospatial Open Source Software GmbH                          #
#                                                                              #
# Contact: Andreas Wytzisk                                                     #
# 52 North Initiative for Geospatial Open Source Software GmbH                 #
# Martin-Luther-King-Weg 24                                                    #
# 48155 Muenster, Germany                                                      #
# info@52north.org                                                             #
#                                                                              #
# This program is free software; you can redistribute and/or modify it under   #
# the terms of the GNU General Public License version 2 as published by the    #
# Free Software Foundation.                                                    #
#                                                                              #
# This program is distributed WITHOUT ANY WARRANTY; even without the implied   #
# WARRANTY OF MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU #
# General Public License for more details.                                     #
#                                                                              #
# You should have received a copy of the GNU General Public License along with #
# this program (see gpl-2.0.txt). If not, write to the Free Software           #
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA or #
# visit the Free Software Foundation web page, http://www.fsf.org.             #
#                                                                              #
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)                          #
# Created: 2010-06-20                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

################################################################################
# PegelOnlineSOS
pegelsos <- SOS(url = "http://v-sos.uni-muenster.de:8080/PegelOnlineSOSv2/sos")
print(object.size(pegelsos), units = c("Mb"))
# works so far... :-)

sosFeaturesOfInterest(pegelsos)

off <- sosOfferings(pegelsos)[[1]]
latestObs <- getObservation(sos = pegelsos,
		offering = off,
#		observedProperty = sosObservedProperties(offering),
		procedure = sosProcedures(off)[11:13],
		latest = TRUE,
		inspect = TRUE) #, verbose = TRUE)
sosResult(latestObs)

# three procedures, but only getting 1 element with one procedure...
pegelsos <- SOS(url = "http://v-sos.uni-muenster.de:8080/PegelOnlineSOSv2/sos")
pegelObs <- getObservation(sos = pegelsos,
		observedProperty = sosObservedProperties(sosOfferings(pegelsos)[[1]])[3],
		offering = sosOfferings(pegelsos)[[1]],
		procedure = sosProcedures(sosOfferings(pegelsos)[[1]])[c(2501,2503,2505)],
		eventTime = sosCreateEventTimeList(time = sosCreateTimePeriod(
						sos = pegelsos,
						begin = Sys.time() - (3600 * 24), # * 360),
						end = Sys.time()))) #, inspect = TRUE)
# Parsing response (size  352 ) ...Finished getObservation to http://v-sos.uni-muenster.de:8080/PegelOnlineSOSv2/sos - received 3 observation(s)/measurement(s) having 87, 2, 2 elements.
# YEAH!

# show parts of the data frame:
pegelObs[[1]]@result[1:5,]

# not enough info? got field descriptions as attributes for each column:
attributes(pegelObs[[1]]@result[,1])
attributes(pegelObs[[1]]@result[,2])
attributes(pegelObs[[1]]@result[,3])


# make plot out of two or three related stations?
range(pegelObs[[1]]@result[,3]); range(pegelObs[[2]]@result[,3])

# Attention: plots ignore the fact that the times do NOT perfectly match!
#x <- 700
#plot(x = obs4[[1]]@result[[1]][1:x], y = obs4[[1]]@result[[3]][1:x], type = "l",
#		col = "steelblue", main = "Temperature in Muenster and Kaernten, 2009",
#		xlab = "Time (00:00 o'clock)",
#		ylab = "Temperature (degree C)",
#		xaxt="n") # do not plot x-axis
#r <- as.POSIXct(round(range(obs4[[1]]@result[[1]]), "days"))
#axis.POSIXct(side = 1, x = obs4[[1]]@result[[1]][1:x], format = "%d. %h",
#		at = seq(r[1], r[2], by="day"))
#lines(x = obs4[[2]]@result[[1]][1:x], y = obs4[[2]]@result[[3]][1:x],
#		col = "orange")
#legend("topleft", legend = c("Muenster", "Kaernten"),
#		col = c("steelblue", "orange"), lty = 1, bty="n")

plot(x = pegelObs[[1]]@result[,1], y = pegelObs[[1]]@result[,3], type = "l")

# Good data?
# Felix's tip: look at the coastal stations, much more interesting!

################################################################################
# Elsterhochwasser, 30.09.2010
procedure <- "Wasserstand-Elster_501390"
offering = sosOfferings(pegelsos)[[2]]
elster <- getObservation(sos = pegelsos, offering = sosOfferings(pegelsos)[[1]],
		procedure = procedure, observedProperty = )
elster[[1]]@result[1:3,]
range(elster[[1]]@result$Wasserstand)
elsterClean <- subset(elster[[1]]@result, Wasserstand > 0)

plot(x = elster[[1]]@result$Time, y = elster[[1]]@result$Wasserstand, ylim = c(100, 600),
		type = "l")

# optional: install the package
#install.packages("sos4R")

# load the sos4R package
#library("sos4R")
#source("/home/daniel/Dropbox/2010_SOS4R/workspace/sos4R/sandbox/loadSources.R")

################################################################################
# pegelonlinesos
# not so nice: not exactly reproducible because data is only stored for 30 days!

# create connection to SOS
pegelsos <- SOS(url = "http://v-sos.uni-muenster.de:8080/PegelOnlineSOSv2/sos")

# what data do I get?
names(sosOfferings(pegelsos))

# set up parameters for request
procedures <- sosProcedures(pegelsos)
procedures <- subset(procedures, procedures %in% grep("*Bake*", procedures, value=TRUE))
# Watch out here because order of elements can change! The correct value is given in the comment.
wasserstand <- sosObservedProperties(pegelsos)[1] # <- Should be "Wasserstand"
wasserstand_roh <- sosOfferings(pegelsos)[[1]] # <- Should be "WASSERSTAND_ROHDATEN"
lastSixtyDays <- sosCreateEventTimeList(time = sosCreateTimePeriod(
				sos = pegelsos,
				begin = Sys.time() - (3600 * 24 * 60),
				end = Sys.time()))

pegelObs <- getObservation(sos = pegelsos, offering = wasserstand_roh,
		observedProperty = wasserstand, procedure = procedures[7],
		eventTime = lastSixtyDays)
data <- sosResult(pegelObs)

# inspect data
summary(data)
data[1:2,]
names(data)
attributes(data[,3])

# clean up data (remove negative values)
data <- subset(data, data[,3]>0)

# create time series from data and plot
library("xts")
bakeA <- xts(x = data[["Wasserstand"]], order.by = data[["Time"]])
plot(bakeA, main = "Water Level at Bake A", xlab = "Time", ylab = "Water Level (cm)", major.ticks = "days")

library("forecast")
# time series plots
tsdisplay(bakeA)
# check the periodicity
periodicity(bakeA)

# fit autoregressive model, selcting complexity by AIC
bakeA.ar <- ar(bakeA)
# forecast and plot
bakeA.ar.fcast <- forecast(bakeA.ar, h = 60 * 48) # 48 hrs
plot(bakeA.ar.fcast)
summary(bakeA.ar.fcast) # overview of forecast
accuracy(bakeA.ar.fcast) # check goodness of fit
# looks interesting! save figure1.png
savePlot(type = "png", filename = "figure1.png")

# TODO add filtered lines

# try out ets
bakeA.ets <- ets(bakeA)
bakeA.ets.fcast <- forecast(bakeA.ets, h = 60 * 24) # 12 hrs
plot(bakeA.ets.fcast)
# not useful

# fit ARMA model by conditional least squares and plot it
bakeA.arma <- arma(bakeA)
plot(bakeA.arma)


################################################################################
# weathersos
#library("sos4R")

# create connection to SOS
weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")

# set up request parameters
stationMuenster <- sosProcedures(weathersos)[[1]]
temperatureOffering <- sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]]
temperature <- sosObservedProperties(weathersos)[5] # "urn:ogc:def:property:OGC::Temperature"
september <- sosCreateTimePeriod(sos = weathersos,
		begin = as.POSIXct("2010-09-01 00:00"),
		end = as.POSIXct("2010-09-30 00:00"))
# make the request
obsSept <- getObservation(sos = weathersos,
		observedProperty = temperature,
		procedure = stationMuenster,
		eventTime = sosCreateEventTimeList(september),
		offering = temperatureOffering)

# inspect data
summary(sosResult(obsSept))
sosResult(obsSept)[1:2,]
names(sosResult(obsSept))
data <- sosResult(obsSept)

# create time series from data and plot
library("xts")
tempSept <- xts(x = data[["urn:ogc:def:property:OGC::Temperature"]],
		order.by = data[["Time"]])
plot(tempSept, main = "Temperature in Muenster",
		xlab = "Time", ylab = "Temperature (degree C)", major.ticks = "weeks")

# time series plots
tsdisplay(tempSept)
# check the periodicity
periodicity(tempSept)
# create models
tempSeptArima <- Arima(tempSept)
tempSeptETS <- ets(tempSept)
# auto-regressive
tempSeptAR <- ar(tempSept)

# create forecasts
library("forecast")
# based on exponential smoothing
fcastETS <- forecast(tempSeptETS, h = 1000)
summary(fcastETS) # overview of forecast
accuracy(fcastETS) # check goodness of fit
plot(fcastETS)
# based on ARIMA model
fcastArima <- forecast(tempSeptArima)
accuracy(fcastArima) # check goodness of fit
plot(fcastArima)


################################################################################
# one year
# create connection to SOS
weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
# set up request parameters
stationMuenster <- sosProcedures(weathersos)[1]
temperatureOffering <- sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]]
temperature <- sosObservedProperties(weathersos)[5] # "urn:ogc:def:property:OGC::Temperature"
timeperiod2009 <- sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
				begin = as.POSIXct("2009-01-01 00:00"),
				end = as.POSIXct("2009-12-31 00:00")))
# make the request
temp2009 <- getObservation(sos = weathersos,
		observedProperty = temperature,
		procedure = stationMuenster,
		eventTime = timeperiod2009,
		offering = temperatureOffering)

# inspect data
summary(sosResult(temp2009))
sosResult(temp2009)[1:2,]
names(sosResult(temp2009))

# create time series from data and plot
library("xts")
tempSeries2009 <- xts(x = sosResult(temp2009)[["urn:ogc:def:property:OGC::Temperature"]],
		order.by = sosResult(temp2009)[["Time"]])
plot(tempSeries2009, main = "Temperature in Muenster", xlab = "Time", ylab = "Temperature (degree C)", major.ticks = "months")

# time series plots
tsdisplay(tempSeries2009)

# check the periodicity
periodicity(tempSeries2009)

# create models
temp2009Arima <- Arima(tempSeries2009)
temp2009ETS <- ets(tempSeries2009)
temp2009AR <- ar(tempSeries2009)
# diagnose model
tsdiag(temp2009ETS)

# create forecasts
library("forecast")

# based on exponential smoothing
fcast2009ETS <- forecast(temp2009ETS, h = 300)
summary(fcast2009ETS) # overview of forecast
accuracy(fcast2009ETS) # check goodness of fit
plot(fcast2009ETS)
plot(fcast2009ETS, xlim=c(26000,30500), xlab=) # last part of the year
# Save figure1.png
savePlot(type = "png", filename = "figure2.png")

# autoregressive model
plot(temp2009AR)

# data is apparently not period enough...
stl(tempSeries2009)
#Fehler in stl(tempSept) : 
#  Zeitreihe ist nicht periodisch oder umfasst weniger als zwei Perioden


################################################################################
# ClimateSOS
climatesos <- SOS("http://giv-sos.uni-muenster.de:8080/ClimateSOS/sos")

length(sosProcedures(climatesos))
# 6

lapply(sosOfferings(climatesos), slot, "name")

################################################################################
# MoodSOS
#
# offerings worth checking out: SummerVillerest, PatientCondition,
# WaterColourNormal, FoamPresence, Microcystin
moodsos <- SOS("http://giv-genesis.uni-muenster.de:8080/52nSOSv3-MoodSOS/sos")

################################################################################
# Umweltbundesamt SOS
umweltsos <- SOS(url = "https://develop.umweltbundesamt.at/SOSsrv/sos")
# https://develop.umweltbundesamt.at/SOSsrv/sos?service=SOS&request=GetCapabilities
# --> 503 Service temporarily unavailable

################################################################################
# Ocean Process Analysis Laboratory, Institute for the Study of Earth, Oceans,
# and Space, University of New Hampshire SOS
COOA_UNH <- SOS("http://www.cooa.unh.edu/cgi-bin/sos/oostethys_sos")
# --> 500 Internal Server Error, and
# --> Capabilities are shown when opening the link above in a browser, but it
# has a strange version: <ows:ServiceTypeVersion>0.0.31</ows:ServiceTypeVersion>

################################################################################
# Gulf of Maine Ocean Observing System SOS
GoMOOS <- SOS("http://www.gomoos.org/cgi-bin/sos/oostethys_sos.cgi",
		version = "0.0.31", verboseOutput = TRUE)
# --> Capabilities are shown when opening the link above in a browser, but it
# has a strange version: <ows:ServiceTypeVersion>0.0.31</ows:ServiceTypeVersion>
# 
# --> Object of class OwsExceptionReport; version: 1.0.0, lang: NA,  1 exceptions (code @ locator : text):
#	MissingParamterValue @ service : No input parameters 

################################################################################
# others:
################################################################################
#
# TAMU <- SOS("http://vastserver.nsstc.uah.edu/vastGC/adcp", verboseOutput = TRUE)
# --> no response
#
# MVCO_WHOI <- SOS("http://mvcodata.whoi.edu/cgi-bin/sos/oostethys_sos")
# --> seems almost empty


################################################################################
# OCEAN STUFF, a lot of interesting data!
#
# http://www.openioos.org/real_time_data/gm_sos.html
#
oceanwatch <- SOS("http://oceanwatch.pfeg.noaa.gov/pysos/sos_mysql2.py",
		method = "GET", verboseOutput = TRUE)
ww6 <- SOS("http://ww6.geoenterpriselab.com:8080/SOS_Weather/sos")

sos-ws <- SOS("http://sos-ws.tamu.edu/tethys/tabs")
# takes forever...		

################################################################################
# some french sos, 52N, but just one week of data....
sandre <- SOS("http://services.sandre.eaufrance.fr/52nSOSv3/sos")

################################################################################
# various
var01.converters <- SosDataFieldConvertingFunctions("urn:terrestris:foss4g:temperature" = sosConvertDouble)
var01 <- SOS(":8280/52nSOSv3_WAR/sos", dataFieldConverters = var01.converters)
time01 <- sosTime(sosOfferings(var01)[[1]])
time01.part <- sosCreateTimePeriod(var01,
		begin = as.POSIXct("2010-08-06"),
		end = as.POSIXct("2010-08-07"))
obs01 <- getObservation(var01, offering = sosOfferings(var01)[[1]],
#		observedProperty = sosObservedProperties(var01),
#		procedure = sosProcedures(var01),
		eventTime = sosCreateEventTimeList(time = time01),
		featureOfInterest = sosCreateFeatureOfInterest(list("8242", "8245")),
		verbose = TRUE)
# feature filter neccessary, otherwise too much data -> works!

var02 <- SOS(":82/cgi-bin/mapserv?map=/tmp/umn/umn_sos.map",
		method = "GET",
		verboseOutput = TRUE)
# no getcapabilities possible:
# 1. url already contains a "?", so if "GET" the error is --- msEvalRegex(): Regular expression error. String failed expression test.
#	-> fixed that case, but then the same error as with "POST"...
# 2. error with "POST" gives HTML page --- Unable to access file. (/tmp/umn/umn_sos.map) 

var03.converters <- SosDataFieldConvertingFunctions(
		"urn:terrestris:foss4g:temperature" = sosConvertDouble,
		"urn:terrestris:foss4g:feature" = sosConvertString)
var03 <- SOS(":8280/deegree-sos-cite100/services",
		dataFieldConverters = var03.converters,
#		method = "POST",
		method = "GET")
#		verboseOutput = TRUE)
var03.off <- sosOfferings(var03)[[1]]
obs03 <- getObservation(sos = var03, offering = var03.off,
		procedure = NA_character_, # must be set, otherwise InvalidParameterValue...
		eventTime = list(NA))
# connection works with GET
# works!
result <- sosResult(obs03)
str(result)
plot(x = result["timestamp"], y = result["val"])
# looks weird because it plots all features

result995 <- subset(x = result, subset = (plz == "995.0\n"),
		select = c(timestamp, val))
plot(x = result995[["timestamp"]], y = result995[["val"]], type = "l")
# works!

# connection with post returns HEX: "3c 3f 78 6d 6c 20 76 65 72 73 "
#  -> translate with http://home2.paulschou.net/tools/xlate/
#  -> ist exception report --- Unable to determine the subcontroller for request type <sos:GetCapabilities...


################################################################################
# SOS at ISE, Verbania, Italy

# add a conversion function for the field definition "...temp"
ise.converters <- SosDataFieldConvertingFunctions("urn:ogc:def:property:OGC:1.0.30:air_temp" = sosConvertDouble)
ise <- SOS("http://sos.ise.cnr.it/sos", dataFieldConverters = ise.converters)
ise
# Offering
ise.offerings <- sosOfferings(ise)

# rain data
ise.rain <- getObservation(sos = ise, offering = ise.offerings[["rain"]], inspect = TRUE)
ise.rain
ise.rain.data <- sosResult(ise.rain)
summary(ise.rain.data)

# getObservation in SOS
ise_air_temp <- getObservation(sos = ise,
		offering = ise.offerings[["air_temperature"]], inspect = TRUE)
ise_air_temp
ise_air_temp.result <- sosResult(ise_air_temp)
summary(ise_air_temp.result)

# request with features of interest
ise.rain.features <- sosFeaturesOfInterest(ise.offerings[["rain"]])

# feature id
ise.rain2 <- getObservation(sos = ise,
		offering = ise.offerings[["rain"]],
		featureOfInterest =  SosFeatureOfInterest(objectIDs = ise.rain.features),
		inspect = TRUE)

# spatial filtering: check out ?sos4R::OGC
bbox <- sosCreateBBOX(lowLat = 10.0, lowLon = 2.0, uppLat = 30.0, uppLon = 5.0,
			srsName = "urn:ogc:def:crs:EPSG:4326")
ise.rain3 <- getObservation(sos = ise,
		offering = ise.offerings[["rain"]],
		featureOfInterest = SosFeatureOfInterest(spatialOps = bbox),
		inspect = TRUE)

################################################################################
# .NET implementation from http://ogc.codeplex.com/ (http://sensordatabus.org/)
# TODO
ws <- SOS("http://ws.sensordatabus.org/Ows/Swe.svc/", method = "GET", verboseOutput = TRUE)
# http://ws.sensordatabus.org/Ows/Swe.svc/?service=SOS&request=GetCapabilities
# Seems not to work with additional request parameters...


################################################################################
# Renaissance Computing Institute
# TODO
renci <- SOS(url = "http://ws.sensordatabus.org/Ows/Swe.svc/", method = "GET")
# does not support sections parameter, capabilities are empty if it's given

# This works:
renci <- SOS(url = "http://ws.sensordatabus.org/Ows/Swe.svc/", method = "GET",
		verboseOutput = TRUE, sections = NA)
renci
renci.off <- sosOfferings(renci)
names(renci.off)

length(renci.off)
# 1772

sosName(renci.off)

################################################################################
# istSOS
# - NO stress tests!
# - Maximum of one month of data!
# - contact: massimiliano.cannata@gmail.com
ist.converters <- SosDataFieldConvertingFunctions(
		"urn:ogc:def:parameter:x-ist::time:iso8601" = sosConvertTime,
		"urn:ogc:def:parameter:x-ist::meteo:air:temperature" = sosConvertDouble)

### here you find 10 minutes aggregated values, they are backward intervals:
# http://geoservice.ist.supsi.ch/sos?request=getcapabilities

# GET #
ist.get <- SOS(url = "http://geoservice.ist.supsi.ch/sos",#verboseOutput = TRUE,
		method = "GET", dataFieldConverters = ist.converters)
sosOfferings(ist.get)

ist.timeperiod <- sosCreateEventTimeList(sosCreateTimePeriod(sos = ist.get,
				begin = as.POSIXct("2011-01-01 00:00"),
				end = as.POSIXct("2011-01-07 00:00")))
ist.offering <- sosOfferings(ist.get)[[1]]

ist.obs <- getObservation(sos = ist.get, verbose = TRUE,
		offering = ist.offering,
#		observedProperty = sosObservedProperties(ist.offering)[2],
		observedProperty = list("urn:ogc:def:parameter:x-ist::meteo:air:temperature"),
		eventTime = ist.timeperiod)
# names of observed properties in capabilities do not match available names:
#Parameter "observedProperty" sent with invalid value: 
#	['urn:ogc:def:property:x-ist::urn:ogc:def:parameter:x-ist::meteo:air:humidity']
# - available options: ['urn:ogc:def:parameter:x-ist::lake:water:height',
#	'urn:ogc:def:parameter:x-ist::meteo:air:humidity',
# [...]

ist.result <- sosResult(ist.obs)
summary(ist.result)
plot(ist.result[["Time"]], ist.result[["air-temperature"]], type = "l")


# POST #
ist.post <- SOS(url = "http://geoservice.ist.supsi.ch/sos",#verboseOutput = TRUE,
		method = "POST", dataFieldConverters = ist.converters)
sosOfferings(ist.post)
ist.obs.post <- getObservation(sos = ist.post, verbose = TRUE,
		offering = ist.offering,
#		observedProperty = sosObservedProperties(ist.offering)[2],
		observedProperty = list("urn:ogc:def:parameter:x-ist::meteo:air:temperature"),
		eventTime = ist.timeperiod)
# ERROR!
# [...] Parameter "offering" is mandatory with multiplicity 1

### raw data:
istraw <- SOS(url = "http://geoservice.ist.supsi.ch/sosraw")
sosOfferings(istraw)

################################################################################
# ENVISION SOS
# TODO
brgm <- SOS(url = "http://swe.brgm.fr/constellation-envision/WS/sos-discovery",
#		verboseOutput = TRUE,
		method = "GET")
summary(brgm)

brgm.offerings <- sosOfferings(brgm)
brgm.offerings
# just one offering

brgm.all <- brgm.offerings[[1]]
sosObservedProperties(brgm.all)
# NONE!

sosProcedures(brgm.all)

################################################################################
# OOSTethys PySOS
# Center for Coastal Margin Observation & Prediction
# http://www.stccmop.org/, http://www.stccmop.org/sos

# TODO implement parsing of om:resultDefinition
# BUT do not do it because the referenced schemas are not available for download
# and the response does not validate with http://schemas.opengis.net/om/1.0.0/observation.xsd

# TODO impmlement method to return data frame
stccmopParseResult <- function(obj, sos, verbose = FALSE) {
	if(verbose) {
		cat("[stccmopParseResult]\n")
		print(obj)
	}
	
	.val <- xmlValue(obj)
	return(.val)
}
print(omResultName)
stccmop <- SOS("http://data.stccmop.org/ws/sos.py", method = "GET",
		parsers = SosParsingFunctions("result" = stccmopParseResult))
sosParsers(stccmop)[[omResultName]]
sosCaps(stccmop)

plot(stccmop) # DOES NOT WORK, missing CRS

names(sosOfferings(stccmop))
sosProcedures(stccmop)

length(sosOfferings(stccmop)); length(unlist(sosProcedures(stccmop)))
# one offering per procedure
unique(unlist(sosObservedProperties(stccmop)))

############################
# Get data for one offering:
off <- sosOfferings(stccmop)[[2]]
sosObservedProperties(off)

# Example from website: http://data.stccmop.org/ws/util/sos.py?service=SOS&request=GetObservation&offering=saturn04&observedProperty=WaterTemperature&responseFormat=text/xml&eventTime=2009-09-26/2009-09-28
getObservation(sos = stccmop, inspect = TRUE, verbose = TRUE,
		observedProperty = list("http://marinemetadata.org/cf#sea_water_temperature"),
		offering = "saturn04")
# not time given returns last observation, fixes in parsing required, see above!

################################################################################
# TODO try out "util" sos @ http://data.stccmop.org/ws/util/sos.py
stccmoputil <- SOS(url = "http://data.stccmop.org/ws/util/sos.py")

################################################################################
# AQE inapplicable error
aqe.converters <- SosDataFieldConvertingFunctions(
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[PM10]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[NO2]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[O3]" = sosConvertDouble,
		"http://www.opengis.net/def/property/OGC/0/SamplingTime" = sosConvertTime,
		"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString)
aqe <- SOS(url = "http://giv-uw.uni-muenster.de:8080/AQE/sos",
		dataFieldConverters = aqe.converters)

prop <- "http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[NO2]"
foi = "foi_DEST080"
time = sosCreateEventTimeList(sosCreateTimePeriod(sos = aqe,
				begin = as.POSIXct("2008-10-01T23:01:00Z"),
		end = as.POSIXct("2008-12-31T23:01:00Z")))
sosOfferings(aqe)[["NO2"]]

obs <- getObservation(sos = aqe, offering = sosOfferings(aqe)[["NO2"]], # verbose = TRUE,
#		observedProperty = sosObservedProperties(sosOfferings(aqe)[["NO2"]]),
		featureOfInterest = sosCreateFeatureOfInterest(list(foi)),
		eventTime = time)
str(obs)
# ist auch bei mir mir inapplicable


################################################################################
# EDC Backup SOS
edc.converters <- SosDataFieldConvertingFunctions(
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[PM10]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[NO2]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[O3]" = sosConvertDouble,
		"http://www.opengis.net/def/property/OGC/0/SamplingTime" = sosConvertTime,
		"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString)

edc <- SOS("http://v-sos.uni-muenster.de:8080/SosAirQuality/sos",
		dataFieldConverters = edc.converters)

sosOfferings(edc)


obs <- getObservation(edc, sosOfferings(edc)[[1]])
result <- sosResult(obs)
summary(result)

################################################################################
# Profiling!
setwd("D:/")
Rprof("EDCprof.out")
obs <- getObservation(edc, sosOfferings(edc)[[1]])
Rprof(NULL)
summaryRprof("EDCprof.out")
# not really useful information

################################################################################
# FH Kaernten
# Contact: Hecke Andreas <A.Hecke@fh-kaernten.at>
cti <- SOS("http://weatherstation.cti.ac.at:8080/52NSOS_CUAS/sos")
cti
summary(cti)
plot(cti)

cti.off <- sosOfferings(cti)
cti.off


cti.off.cuas <- sosOfferings(cti)[["WEATHER_CUAS"]]

time = sosCreateEventTimeList(sosCreateTimePeriod(sos = cti,
				begin = as.POSIXct(Sys.time() - 3600),
				end = as.POSIXct(Sys.time())))

cuas.lastHour <- getObservation(cti, cti.off.cuas, verbose = TRUE,
#		inspect = TRUE,
		eventTime = time)
#[.sosRequest_1.0.0] response:
#raw(0)
#raw as char:   
#Error in if (regexpr("(<html>|<HTML>|<!DOCTYPE HTML)", .response) > 0) { : 
#  argument is of length zero

# Same request pasted in to basic test client gives result!

# Try with get... even weirder.
ctiget <- SOS("http://weatherstation.cti.ac.at:8080/52NSOS_CUAS/sos",
		method = "GET")
getObservation(ctiget, cti.off.cuas, verbose = TRUE,
#		inspect = TRUE,
		eventTime = time)
#Object of class OwsExceptionReport; version: 1.0.0; lang: NA;
#1 exception(s) (code @ locator : text):
#		InvalidRequest @ REQUEST :
#		The GET request GetObservation is not supported by this SOS. 

#
#
#
cuas.lastHour <- getObservation(cti, cti.off.cuas, verbose = TRUE,
#		inspect = TRUE,
		procedure = sosProcedures(cti.off.cuas)[[1]],
		eventTime = time)

################################################################################
# some services from the AGILE 2011 paper "Empirical Study ..."
# TODO test, create demos if interesting data
# use time series analysis like presented at
# http://www.r-bloggers.com/time-series-analysis-and-mining-with-r/

# WEATHERFLOW
weatherflow <- SOS(url = "http://www.weatherflow.com/sos/sos.pl")

# offering network-all
# TODO make demo with weatherflow

# WAVCIS (part of MMI)
# http://www.wavcis.lsu.edu/SOS/server.asp?request=GetCapabilities
wavcis <- SOS(url = "http://www.wavcis.lsu.edu/SOS/server.asp")


# Sensor Data Bus 
# http://www.sensordatabus.org/Pages/SOS.aspx
# http://ogc.codeplex.com/
# http://ws.sensordatabus.org/Ows/Swe.svc/?service=SOS&request=GetCapabilities
sdb <- SOS(url = "http://ws.sensordatabus.org/Ows/Swe.svc/")

# Pegelonline @ WSV
pegelonline <- SOS(url = "http://www.pegelonline.wsv.de/webservices/gis/sos")

# ccip projct, degree SOS provided by Lat-Lon
# http://ccip.lat-lon.de/
# http://ccip.lat-lon.de/ccip-sos/services?request=GetCapabilities&service=SOS
ccip <- SOS(url = "http://ccip.lat-lon.de/ccip-sos/services")

# CCIW
# http://devgeo.cciw.ca/cgi-bin/mapserv/sostest?request=GetCapabilities?service=SOS
# Supposed to be Mapserv SOS, but strange error messages about WMS and WFS...

# CSE, probably 52N
# http://sensorweb.cse.unt.edu:8080/teo/sos?request=GetCapabilities&service=SOS

# DISL
# Dolphins!
disl <- SOS(url = "http://gcoos.disl.org/cgi-bin/oostethys_sos.cgi")

# RSMAS
# gulf of mexico!
# http://gcoos.rsmas.miami.edu/sos_server.php?service=SOS&request=GetCapabilities
rsmas <- SOS(url = "http://gcoos.rsmas.miami.edu/sos_server.php")

# INESCPORTO
# several example services for oostethys and pysos
# http://gis.inescporto.pt/

# HIOOS
# Hawaii!
# http://oos.soest.hawaii.edu/oostethys/sos?service=SOS&request=GetCapabilities
hioos <- SOS(url = "http://oos.soest.hawaii.edu/oostethys/sos")

# NASA Real Time Mission Monitor
# TODO real interesting data, check it out!
# http://rtmm.nsstc.nasa.gov/
# UAH VAST SOS
# http://rtmm2.nsstc.nasa.gov/SOS/footprint?request=GetCapabilities&service=SOS&version=1.0.0
footprint <- SOS(url = "http://rtmm2.nsstc.nasa.gov/SOS/footprint")
# http://rtmm2.nsstc.nasa.gov/SOS/nadir?request=GetCapabilities&service=SOS&version=1.0.0
nadir <- SOS(url = "http://rtmm2.nsstc.nasa.gov/SOS/nadir")


################################################################################
# http://sensorweb.forum.52north.org/52-North-thin-client-and-istSOS-tp3370604p3370604.html
# TODO check it out!
#
sadco.url <- "http://ict4eo.meraka.csir.co.za/sadcosos/sos.py"
sadco <- SOS(url = sadco.url, verboseOutput = TRUE)
# problem: missing elements in one ObservationOffering so the object cannot be
# created:
#Error in validObject(.Object) : 
#		invalid class "SosObservationOffering" object: invalid object for slot "procedure" in class "SosObservationOffering": got class "list", should be or extend class "character"
#In addition: Warning messages:
#		1: In FUN(X[[1L]], ...) :
#		Mandatory element 'responseFormat' missing in offering temporary
#2: In FUN(X[[1L]], ...) :
#		Mandatory element 'responseMode' missing in offering temporary
#3: In FUN(X[[1L]], ...) :
#		Mandatory element 'time' missing in offeringtemporary
xmlstring <- '<sos:ObservationOffering gml:id="temporary">
		<gml:name>urn:x-sadco::offering:temporary</gml:name>
		<gml:description>temporary offering to hold self-registered procedures/sensors waiting for service adimistration acceptance</gml:description>
		<gml:boundedBy>
		<gml:null>inapplicable</gml:null>
		</gml:boundedBy>
		</sos:ObservationOffering>'
offering <- xmlParseString(xmlstring)

weathersos <- SOS(SosExampleServices()[[1]])
weathersos@verboseOutput <- TRUE

.procedure <- offering[sosProcedureName]
if(length(.procedure) < 1) {
	print("lala")
	.procedure <- as.character(c())
}
str(.procedure)

offering_parsed <- parseSosObservationOffering(offering, weathersos)
str(offering_parsed)
# fixed it!

sadco <- SOS(url = sadco.url, verboseOutput = TRUE)
# WORKS with a few warning messages

sosCapabilitiesDocumentOriginal(sos=sadco)
sosContents(sadco)

getObservation(sadco, offering = sosOfferings(sadco)[[2]],
		eventTime = sosCreateTime(sadco, "2010-01-01::2010-01-05"))
# must be some error in SOS... response:
#
#Object of class OwsExceptionReport; version: 1.0.0; lang: NA;
#1 exception(s) (code @ locator : text):
#		1 @ service :
#		
#		Parameter "offering" is mandatory with multiplicity 1
#
# but there is an offering in the request!
# TODO check getting data from service again

################################################################################
# Spatial Data Infrastructure of the CNR-Institute for Atmospheric Pollution
# - http://sdi.iia.cnr.it/geoint/
# - is a 52N SOS
library("sos4R")
iia <- SOS("http://sdi.iia.cnr.it/sos/sos")

# has 52N dummy data
sosOfferings(iia)

# check out offering "SEA CRUISE 2003"
seacruise <- sosOfferings(iia)["SEA CRUISE 2003"]

# TODO plot it, use spacetime classes


################################################################################
# TODO check out SOSs in ERDDAP
# these are NOAA SOSs, but maybe new ones, you never know...
#
# http://coastwatch.pfeg.noaa.gov/erddap/info/index.html

################################################################################
# TODO check out Wupperverband SOS
#
library("sos4R")
fluggs <- SOS(url = "http://fluggs.wupperverband.de/sos/sos")

sosContents(fluggs)

sosProcedures(fluggs)
sosObservedProperties(fluggs)


################################################################################
# TODO try out new CO-OPS SOSs
# 
#CO-OPS has expanded its SOS services with an addition of the following 7 new services.
#
#One Minute Water Level Data
#Six Minute Water Level Data
#Hourly Height Water Level Data
#High Low Water Level Data
#Daily Mean Water Level Data
#Harmonic Constituents
#Datums

#These services are offered for single station and as collections.
#In addition, CO-OPS is now providing its observational data in KML format. All
#CO-OPS SOS services, including the newly added ones, can be retrieved in KML.
#
#Please note that these service are presently available on the evaluation test
#site (http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos-test/) till November
#30th 2011.
#
#On December 1, 2011 at 10 am EDT, CO-OPS will add these new changes to our
#operational SOS web site (http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos/).
#A reminder email will be sent out that week.

ioos_testing <- SOS(url = "http://opendap.co-ops.nos.noaa.gov/ioos-dif-sos-test/SOS",
#		method = SosSupportedConnectionMethods()[["GET"]]
)
sosObservedProperties(ioos_testing)[["network-All"]]
unique(unlist(sosObservedProperties(ioos_testing)))

################################################################################
# TODO check out SOS from Sandre, French National Service for Water Data and 
#                         Common Repositories Management 
#                         http://sandre.eaufrance.fr/
#
# WaterML response format!
#
# http://services.sandre.eaufrance.fr/52nSOSv3_WML/sos?request=GetCapabilities&service=SOS
library("sos4R")

sandre_converters <- SosDataFieldConvertingFunctions(
	"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString)
sandre <- SOS(url = "http://services.sandre.eaufrance.fr/52nSOSv3_WML/sos",
		dataFieldConverters = sandre_converters)

sosObservedProperties(sandre)
sosProcedures(sandre)

sosResponseFormats(sandre)$GetObservation
# not mentioning WaterML, probably OK
sosResultModels(sandre)[["GetObservation"]]
# wml2:TimeseriesObservation !
#
# TODO try to request WaterML TimeseriesObservation, must add the namespace to
# the request and need a mechanism for that.
# TODO parse WaterML TimeseriesObservation
#

##########
# GET DATA
myOffering <- sosOfferings(sandre)[["A2140100"]]

# TODO the time handling workflow based on the available time should be easier
# TODO implement sosTime with convert option, possible default it to TRUE?
lastTime <- sosTime(sosTime(myOffering)@endPosition)
lastTime.posix <- as.POSIXct(lastTime)
startTime.posix <- lastTime.posix - 3600 * 24 * 7 # one week
myTimeString <- paste(startTime.posix, "::", lastTime.posix, sep = "")
myTime <- sosCreateTime(sos = sandre, time = myTimeString)
myTime

lastDayA2140100 <- getObservation(sos = sandre, offering = myOffering,
		eventTime = myTime)
summary(sosResult(lastDayA2140100[[1]]))
summary(sosResult(lastDayA2140100[[2]]))

sosResult(lastDayA2140100)
# cannot join because there are two differen observations

waterheight <- sosResult(lastDayA2140100[[2]])
plot(waterheight)

plot(x = waterheight$SamplingTime, y = waterheight$WATERHEIGHT, type = "l",
		main = paste(sosId(myOffering), sosName(myOffering)),
		# TODO should work with sosName(waterheight) or better sosName(lastDayA2140100[[1]])
		ylab = paste("Waterheight [", sosUOM(waterheight), "]"),
		xlab = "sampling time")



############
# PROCEDURES
# testing handling of multiple sensors in describeSensor
sensor_1_1 <- describeSensor(sos = sandre, # verbose = TRUE,
		procedure = sosProcedures(sandre)[[1]])
str(sensor_1_1)
# not discovery profile, no elements are found:
sensor_1_1[[1]]@xml

# multiple sensors work, but saving, too?
getwd()
describeSensor(sos = sandre, procedure = sosProcedures(sandre)[[1]],
		saveOriginal = TRUE)

sosCoordinates(sensor_1_1)
plot(sensor_1_1) # does not work... yet

#
# DATA
#
sosOfferings(sandre)[[1]] # check valid time interval
myTime <- sosCreateTime(sos = sandre, time = "2011-10-18::2011-10-20")
myOffering <- sosOfferings(sandre)[[1]]
obs_1 <- getObservation(sos = sandre, verbose = TRUE,
		offering = myOffering, inspect = TRUE,
		procedure = sosProcedures(myOffering)[[1]],
		# limit to one procedure during testing, works
		eventTime = myTime)
# first time run with warnings for missing converters, added to defaults, was
# http://www.opengis.net/def/property/OGC/0/SamplingTime and
str(obs_1)

#################
#Warning message:
#		In if (.contentType == mimeTypeXML) { :
#					the condition has length > 1 and only the first element will be used


# TODO fix observed properties for OmObservation and OmObservationCollection
sosObservedProperties(obs_1_[[1]])
obs_1[[1]]@observedProperty


obs <- getObservation(sos = sandre, verbose = TRUE,
		offering = myOffering, inspect = TRUE,
		eventTime = myTime)
#######
# Error in match.names(clabs, names(xi)) : 
#  names do not match previous names
# FIXME

################################################################################
# TODO check out SOS from Spain with environmental data
# http://elcano.dlsi.uji.es:8082/SOSM2/sos?request=GetCapabilities&service=SOS
elcano <- SOS(url = "http://elcano.dlsi.uji.es:8082/SOSM2/sos")

sosContents(elcano)

sosOfferings(elcano)

################################################################################
# TODO check out SOS from Sensors4All Framework, FH Kärnten
sensors4all <- SOS(url = "http://weatherstation.cti.ac.at:8080/52nSOSv3_CUAS/sos")

sosOfferings(sensors4all)
sensors4all_off <- sosOfferings(sensors4all)[[1]]
plot(sensors4all)

library(maps); library(mapdata); library(maptools)
if(!require(rgdal, quietly = TRUE))
	print("rgdal not present: CRS values will not be converted correctly")
data(worldHiresMapEnv)
crs <- sosGetCRS(sensors4all)
region <- map.where(database = "worldHires",
		sosCoordinates(sensors4all_off)) # find region
worldHigh <- pruneMap(map(database = "worldHires", region = region,
				plot = FALSE))
worldHigh_Lines <- map2SpatialLines(worldHigh, proj4string = crs)
plot(worldHigh_Lines, col = "grey50")
plot(sensors4all_off, add = TRUE, lwd = 3)


################################################################################
# TODO check out VITO SOS with air quality data (live data planned) for Antwerp
library(sos4R)
idea <- SOS(url = "http://sensorweb.vito.be:8080/IDEA_52nSOSv3.2/sos")

idea_offs <- sosOfferings(idea)
sosObservedProperties(idea)

alix5_off_time <- sosTime(idea_offs$alix5)
alix5_time <- sosCreateEventTimeList(alix5_off_time)
#myTime <- sosCreateTime(sos = idea, time = "2011-10-18::2011-10-20")

alix5_observation <- getObservation(sos = idea, offering = idea_offs$alix5,
		eventTime = alix5_time, verbose = FALSE)
alix5_result <- sosResult(alix5_observation)


################################################################################
# TODO check out IRCEL - CELINE SOS with air quality data
library(sos4R)
ircel <- SOS(url = "http://sos.irceline.be/sos")

sosOfferings(ircel)
