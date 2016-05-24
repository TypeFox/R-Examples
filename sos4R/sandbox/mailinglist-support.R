################################################################################
# This script file collects questions from the mailing list and forum          #
# available at http://geostatistics.forum.52north.org/.                        #
#                                                                              #
# Please note that code is supplied by the respective mailing list or forum    #
# member and only copied here for testing.                                     #
################################################################################

################################################################################
# Re: [52N Geostatistics] about sos4R
# Mon, 22 Nov 2010 10:37:34 +0100 (CET)
# http://geostatistics.forum.52north.org/Re-52N-Geostatistics-about-sos4R-tp1909723p1944283.html

# add a conversion function for the field definition "...:chla_conc"
ise_chla.converters <-
		SosDataFieldConvertingFunctions("urn:ogc:def:property:OGC:1.0.30:chla_conc"
						= sosConvertDouble)
ise <- SOS("http://sos.ise.cnr.it/sos", dataFieldConverters =
				ise_chla.converters)
ise

# Offering
ise.offerings <- sosOfferings(ise)
ise.offerings

# set up request parameters
sosProcedures(ise)
station_ise_chla <- sosProcedures(ise)[[5]] # [[5]] mean second procedures:
total_chla_FP
station_ise_chla
sosOfferings(ise)
Offering_ise_chla <- sosOfferings(ise)[["total_chla_FP"]]
Offering_ise_chla
ise.chla.features <- sosFeaturesOfInterest(ise.offerings[["total_chla_FP"]])
# spatial filtering: check out 
ise.chla.features <- SosFeatureOfInterest(objectIDs =
				ise.chla.features[1:11])
ise.chla.features
ise_chla_time <- sosCreateTimePeriod(sos = ise,
		begin = as.POSIXct("2005-08-30 00:00"),
		end = as.POSIXct("2005-08-30 23:00"))
ise_chla_time

# make the request
obs_ise_chla <- getObservation(sos = ise,
		procedure = station_ise_chla,
		offering = Offering_ise_chla,
		featureOfInterest = ise.chla.features,
		eventTime = sosCreateEventTimeList(ise_chla_time),
		inspect = TRUE)

data_ise_chla <- sosResult(obs_ise_chla)
data_ise_chla

# extract the positions
coords_ise_chla <- sosCoordinates(obs_ise_chla)

# attach the positions to the data frame (attention: matching names different in SVN version)
data_coords_ise_chla <- merge(x = data_ise_chla, y = coords_ise_chla,
		by.x = "feature", by.y = "foi_id")

# after version 0.1-07:
data_coords_ise_chla <- merge(x = data_ise_chla, y = coords_ise_chla)

# some example rows of the merged data.frame
data_coords_ise_chla[c(1:3, 101:103, 301:303),]


#
# 2011-07-13: Support for a direct question
#
# 1. I managed to connect to the following server, but when I tried the spanish one, I get the error below
#sos <- SOS(url="http://elcano.dlsi.uji.es:8080/SOS_MCLIMATIC")
#Error in from:to : result would be too long a vector
#In addition: Warning messages:
#		1: In strsplit(str, "\\\r\\\n") : input string 1 is invalid in this locale
#2: In grep("^[[:space:]]+$", str) :
#		input string 1 is invalid in this locale
#3: In grep("^HTTP", lines) : input string 1 is invalid in this locale
#4: In max(i) : no non-missing arguments to max; returning -Inf
#5: In parseHTTPHeader(c(header, str)) : NAs introduced by coercion
#6: In max(i) : no non-missing arguments to max; returning -Inf

mclimatic <- SOS(url="http://elcano.dlsi.uji.es:8080/SOS_MCLIMATIC/sos")
# SOLUTION: wrong endpoint, missing /sos
summary(mclimatic)



################################################################################
# [52N Geostatistics] SOS error
# 01.12.2011 17:50
# http://geostatistics.forum.52north.org/SOS-error-td3552095.html

library(sos4R)
library(sp)
library(rgdal)
watershapSOSConverters <-
		SosDataFieldConvertingFunctions(
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22045DS"=
						sosConvertDouble,
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22047DS"=
						sosConvertDouble,
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22020DS"=
						sosConvertDouble,
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22020US"=
						sosConvertDouble,
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22045US"=
						sosConvertDouble,
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22046DS"=
						sosConvertDouble,
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22019DS"=
						sosConvertDouble,
				"urn:ogc:object:feature:sensor:VE:waterlevel-sensor-22046US"=
						sosConvertDouble,
				"http://www.opengis.net/def/property/OGC/0/SamplingTime" = 
						sosConvertTime,
				"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" =
						sosConvertString)
watershapSOS <-
		SOS(url="http://137.224.18.29:8080/WS_VenE_SOSv3/sos",
				dataFieldConverters = watershapSOSConverters)
summary(watershapSOS)
watershapOffering <- sosOfferings(watershapSOS)
watershapTimeSOS= sosCreateTime(sos = watershapSOS, time =
				"2011-10-03::2011-11-10")
waterLevel<-watershapOffering[["WATERLEVEL"]]
obsWaterLevel=getObservation(sos=watershapSOS, offering=waterLevel,
		eventTime=watershapTimeSOS)
WaterLevelGauges=sosResult(obsWaterLevel, coordinates = TRUE)
##Get CRS of Observations
obsWaterLevelCRS<-sosGetCRS(obsWaterLevel)
##Create spatial objects from observations
coordinates(WaterLevelGauges)=~lat+lon
proj4string(WaterLevelGauges)=obsWaterLevelCRS
##Set Coordinate System Amersfoort
Amersfoort.RD.New=CRS("+proj=sterea +lat_0=52.15616055555555
				+lon_0=5.38763888888889 +k=0.999908 +x_0=155000 +y_0=463000 +ellps=bessel
				+units=m
				+towgs84=565.2369,50.0087,465.658,-0.406857330322398,0.350732676542563,-1.8703473836068,4.0812
				+no_defs +no_defs")
WaterLevelGaugesProj<-spTransform(WaterLevelGauges,Amersfoort.RD.New)

#*Everything seems to be OK, but I am receiving an error message with the
#following text*
#		
#[sos4R] Received response (size: 418928 bytes), parsing ...
#[sos4R] Finished getObservation to
#http://137.224.18.29:8080/WS_VenE_SOSv3/sos 
#--> received 7 observation(s) having 9041 result values [ 1380, 1138, 1138,
#		1300, 1301, 1392, 1392 ]. 
#Warning message:
#		In if (.contentType == mimeTypeXML) { :
#					the condition has length > 1 and only the first element will be used

# DN:
plot(WaterLevelGaugesProj)

getObservation(sos = watershapSOS, offering = waterLevel,
		# inspect = TRUE,
		verbose = TRUE,
		eventTime = watershapTimeSOS)

# SOLUTION:
# Error message can be ignored, warning is handled with verbose message in 0.2-7


################################################################################
# [Bug 697]
# https://bugzilla.52north.org/show_bug.cgi?id=697
#
# Summary: [sos4R : getObservation] Erreur dans sum(.resultLength) : 'type'
#		(list) de l'argument incorrect
#		Product: 52N Geostatistics
#		Version: unspecified
#		Platform: PC
#		OS/Version: Windows
#		Status: NEW
#		Severity: enhancement
#		Priority: Lowest
#		Component: sos4R
#		AssignedTo: d.nuest@52north.org
#		ReportedBy: guillaume.wattelez@univ-nc.nc
#		Estimated Hours: 0.0

#:sos4R package:
#		Overview : Error with the function getObservation of the R package sos4R on my
#SOS DataBase and the eo2h DataBase.

#Steps to reproduce : In R.
library(sos4R);
eo2h = SOS(url = "http://141.30.100.135:8080/eo2heavenSOS/sos");

eo2h_off = sosOfferings(eo2h);
obs = getObservation(eo2h, offering = eo2h_off[[16]]);
obs = getObservation(eo2h, offering = eo2h_off[[14]]);
#Erreur dans sum(.resultLength) : 'type' (list) de l'argument incorrect
#		De plus : Il y a eu 50 avis ou plus (utilisez warnings() pour voir les 50
#		premiers)
#		>     # There's an error #

#Description : 
#		I have the same error when I try with my SOS DataBase
#Actual Result : 
#		"[sos4R] Received response (size: 2746928 bytes), parsing ...
#		Erreur dans sum(.resultLength) : 'type' (list) de l'argument incorrect"
#with some warnings (put at the end).
#Finally, I don't know if the problem is from R or the SOS DataBase filling. If
#		it can help, a xml response is send when I request a getObservation in my
#		browser.

# FIX:

# First, fix the warnings, see demo("eo2heaven")
eo2h_converters <- SosDataFieldConvertingFunctions(
		"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString,
		"http://www.eo2heaven.org/classifier/parameter/daily_average/BEN" = sosConvertDouble)
eo2h = SOS(url = "http://141.30.100.135:8080/eo2heavenSOS/sos",
		dataFieldConverter = eo2h_converters);

eo2h_off = sosOfferings(eo2h);
obs16 = getObservation(eo2h, offering = eo2h_off[[16]]); # NOX
# works fine, no warnings
summary(obs16)
obs16_data <- sosResult(obs16)
summary(obs16_data)
# all OK!

eo2h_off[[14]] # benzene
obs14 = getObservation(eo2h, offering = eo2h_off[[14]]);
# again, warnings > add another field to converters
# actually fixed in next release because of adding ug/m3 to known unit list.

# Solution: None, the error does not occur on my system!
