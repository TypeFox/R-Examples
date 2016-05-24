### R code from vignette source 'sos4R.Rnw'

###################################################
### code chunk number 1: options
###################################################
.goOnline <- FALSE # triggers whether data is downloaded and stored in /vignettes folder
.verbose <- TRUE

#library("cacheSweave")
#setCacheDir("/tmp/cacheSweave/sos4R")

.getFilePath <- function(name) {
	if(regexpr(pattern = ".RData", text = name) > 0)
		.fileEnding <- ""
	else .fileEnding <- ".xml"

	# since Vignette is in /vignettes/ this doesn't seem to work anymore...
#	.path <- paste(find.package("sos4R", lib.loc = NULL), "/vignettes/",
#			name, .fileEnding, sep = "")
	.path <- paste0(name, .fileEnding)
	if(.verbose) cat("Loading file", .path, "\n")
	
	return(.path)
}


###################################################
### code chunk number 2: sos4R.Rnw:65-68
###################################################
options(width=60)
options(SweaveHooks=list(fig=function()
					par(mar=c(5.1, 4.1, 1.1, 2.1))))


###################################################
### code chunk number 3: load
###################################################
library("sos4R")


###################################################
### code chunk number 4: supported01
###################################################
SosSupportedOperations()
SosSupportedServiceVersions()
SosSupportedConnectionMethods()
SosSupportedResponseFormats()


###################################################
### code chunk number 5: supported02
###################################################
SosSupportedResponseModes()
SosSupportedResultModels()


###################################################
### code chunk number 6: supported03 (eval = FALSE)
###################################################
## SosSupportedSpatialOperators()


###################################################
### code chunk number 7: supported04
###################################################
toString(SosSupportedSpatialOperators())


###################################################
### code chunk number 8: supported05 (eval = FALSE)
###################################################
## SosSupportedTemporalOperators()


###################################################
### code chunk number 9: supported06
###################################################
toString(SosSupportedTemporalOperators())


###################################################
### code chunk number 10: default
###################################################
SosDefaultConnectionMethod()
SosDefaults()


###################################################
### code chunk number 11: converterFunc (eval = FALSE)
###################################################
## SosEncodingFunctions()
## SosParsingFunctions()
## SosDataFieldConvertingFunctions()


###################################################
### code chunk number 12: conn
###################################################
.fileMySOS <- "mySOS.RData" 
if(.goOnline) {
	mySOS <- SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
	save(mySOS, file = .fileMySOS)
} else {
	load(.fileMySOS)
}


###################################################
### code chunk number 13: conn (eval = FALSE)
###################################################
## mySOS <- SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")


###################################################
### code chunk number 14: connDetails1 (eval = FALSE)
###################################################
## sosUrl(mySOS)
## sosTitle(mySOS)
## sosAbstract(mySOS)
## sosVersion(mySOS)
## sosTimeFormat(mySOS)
## sosMethod(mySOS)


###################################################
### code chunk number 15: connDetails2 (eval = FALSE)
###################################################
## sosEncoders(mySOS)
## sosParsers(mySOS)
## sosDataFieldConverters(mySOS)


###################################################
### code chunk number 16: connDetails3
###################################################
mySOS
summary(mySOS)


###################################################
### code chunk number 17: capsOriginal (eval = FALSE)
###################################################
## sosCapabilitiesDocumentOriginal(sos = mySOS)


###################################################
### code chunk number 18: getCap1 (eval = FALSE)
###################################################
## getCapabilities(sos = mySOS)


###################################################
### code chunk number 19: getCap2 (eval = FALSE)
###################################################
## sosServiceIdentification(mySOS)
## sosServiceProvider(mySOS)
## sosFilter_Capabilities(mySOS)
## sosContents(mySOS)


###################################################
### code chunk number 20: getCap3
###################################################
sosTime(mySOS)


###################################################
### code chunk number 21: getCap4 (eval = FALSE)
###################################################
## sosOperationsMetadata(mySOS)
## sosOperation(mySOS, "GetCapabilities")
## sosOperation(mySOS, sosGetCapabilitiesName)


###################################################
### code chunk number 22: getCap5 (eval = FALSE)
###################################################
## sosResponseFormats(mySOS)
## sosResponseMode(mySOS)
## sosResultModels(mySOS)


###################################################
### code chunk number 23: getCap6a (eval = FALSE)
###################################################
## sosResponseMode(mySOS, unique = TRUE)


###################################################
### code chunk number 24: getCap6b
###################################################
toString(sosResponseMode(mySOS, unique = TRUE))


###################################################
### code chunk number 25: getCap7
###################################################
sosResultModels(mySOS)[1:3]


###################################################
### code chunk number 26: getCap8a (eval = FALSE)
###################################################
## sosResponseMode(mySOS)[[sosGetObservationByIdName]]


###################################################
### code chunk number 27: getCap8b
###################################################
toString(sosResponseMode(mySOS)[[sosGetObservationByIdName]])


###################################################
### code chunk number 28: getCap9a (eval = FALSE)
###################################################
## sosResultModels(mySOS)[[sosGetObservationName]][3:4]


###################################################
### code chunk number 29: getCap9b
###################################################
toString(sosResultModels(mySOS)[[sosGetObservationName]])


###################################################
### code chunk number 30: getCap10 (eval = FALSE)
###################################################
## sosResponseFormats(mySOS)[[sosGetObservationByIdName]]


###################################################
### code chunk number 31: getCap10
###################################################
toString(paste(sosResponseFormats(mySOS)[[sosGetObservationByIdName]]))


###################################################
### code chunk number 32: sosGetCRS
###################################################
sosGetCRS("urn:ogc:def:crs:EPSG:4326")

# returns the CRS of offering(s) based on the CRS 
# used in the element gml:boundedBy:
sosGetCRS(mySOS)[1:2]

sosGetCRS(sosOfferings(mySOS)[[1]])


###################################################
### code chunk number 33: sosPlot
###################################################
# background map:
library(maps); library(mapdata); library(maptools)
data(worldHiresMapEnv)
crs <- sosGetCRS(mySOS)[[1]]
worldHigh <- pruneMap(
		map(database = "worldHires",
			region = c("Germany", "Austria", "Netherlands"),
			plot = FALSE))
worldHigh.lines <- map2SpatialLines(worldHigh, proj4string = crs)

# the plot:
plot(worldHigh.lines, col = "grey50")
plot(mySOS, add = TRUE, lwd = 3)
title(main = paste("Offerings by '", sosTitle(mySOS), "'", sep = ""),
		sub = toString(names(sosOfferings(mySOS))))


###################################################
### code chunk number 34: sosPlotFigure
###################################################
getOption("SweaveHooks")[["fig"]]()
# background map:
library(maps); library(mapdata); library(maptools)
data(worldHiresMapEnv)
crs <- sosGetCRS(mySOS)[[1]]
worldHigh <- pruneMap(
		map(database = "worldHires",
			region = c("Germany", "Austria", "Netherlands"),
			plot = FALSE))
worldHigh.lines <- map2SpatialLines(worldHigh, proj4string = crs)

# the plot:
plot(worldHigh.lines, col = "grey50")
plot(mySOS, add = TRUE, lwd = 3)
title(main = paste("Offerings by '", sosTitle(mySOS), "'", sep = ""),
		sub = toString(names(sosOfferings(mySOS))))


###################################################
### code chunk number 35: defaultValue
###################################################
defaultOutputFormatDescribeSensor <- gsub(pattern = "&quot;", replacement = "'", x = sosDefaultDescribeSensorOutputFormat)


###################################################
### code chunk number 36: describeSensor1a
###################################################
.sensorFile <- "mySensor"
if (.goOnline) {
	mySensor <- describeSensor(sos = mySOS, # verbose = TRUE
			procedure = "urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93",
			saveOriginal = .sensorFile)
} else {
	# mySensor <- parseFile(mySOS, .getFilePath(.sensorFile), verbose = TRUE)
	#mySensor <- parseFile(mySOS, paste0(.sensorFile, ".xml")) #, verbose = TRUE)
	mySensor <- parseFile(mySOS, .getFilePath(.sensorFile), verbose = .verbose)
}


###################################################
### code chunk number 37: describeSensor1b (eval = FALSE)
###################################################
## # manual assignment used because procedure order might change:
## mySensor <- describeSensor(sos = mySOS,
## 		procedure = "urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93")
## 
## # using procedure referencing:
## mySensor <- describeSensor(sos = mySOS,
## 		procedure = sosProcedures(obj = mySOS)[[1]][[1]])


###################################################
### code chunk number 38: describeSensor1c
###################################################
mySensor


###################################################
### code chunk number 39: describeSensor2 (eval = FALSE)
###################################################
## sosCoordinates(mySensor)


###################################################
### code chunk number 40: describeSensor3
###################################################
sosId(mySensor)
sosName(mySensor)
sosAbstract(mySensor)


###################################################
### code chunk number 41: describeSensorCoords
###################################################
mySensor.coords <- sosCoordinates(mySensor)
attributes(mySensor.coords)
sosBoundedBy(mySensor)


###################################################
### code chunk number 42: describeSensorPlotCodeForText (eval = FALSE)
###################################################
## library(maps); library(mapdata); library(maptools)
## data(worldHiresMapEnv)
## 
## # get sensor descriptions
## procs <- unique(unlist(sosProcedures(mySOS)))
## procs.descr <- lapply(X = procs, FUN = describeSensor, sos = mySOS)
## 
## sensors.crs <- unique(sosGetCRS(procs.descr))[[1]]
## worldHigh <- pruneMap(map(database = "worldHires",
## 				region = c("Germany", "Austria", "Netherlands",
## 						"Italy"),
## 				plot = FALSE))
## worldHigh.lines <- map2SpatialLines(worldHigh, proj4string = sensors.crs)
## 
## plot(worldHigh.lines, col = "grey50", ylim = c(44.0, 54.8))
## for(x in procs.descr)
## 	plot(x, add = TRUE, pch = 19)
## text(sosCoordinates(procs.descr)[c("x", "y")],
## 		labels = sosId(procs.descr), pos = 4, cex = 0.8)
## title(main = paste("Sensors of", sosTitle(mySOS)))


###################################################
### code chunk number 43: describeSensorPlot
###################################################
getOption("SweaveHooks")[["fig"]]()
library(maps); library(mapdata); library(maptools)
data(worldHiresMapEnv)

# get sensor descriptions
.fileProcs <- "procs.descr.RData"
if (.goOnline) {
	procs <- unique(unlist(sosProcedures(mySOS)))
	procs.descr <- lapply(X = procs, FUN = describeSensor, sos = mySOS)
	save(procs.descr, file = .fileProcs)
} else {
	load(.fileProcs)
}

sensors.crs <- unique(sosGetCRS(procs.descr))[[1]]
worldHigh <- pruneMap(map(database = "worldHires",
				region = c("Germany", "Austria", "Netherlands",
						"Italy"),
				plot = FALSE))
worldHigh.lines <- map2SpatialLines(worldHigh, proj4string = sensors.crs)

plot(worldHigh.lines, col = "grey50", ylim = c(44, 55))
# only plot the procedures where I can find the id (and then probably the
# position can also be found
for(x in procs.descr)
	plot(x, add = TRUE, pch = 19)
text(sosCoordinates(procs.descr)[c("x", "y")],
		labels = sosId(procs.descr), pos = 4, cex = 0.6)
title(main = paste("Sensors of", sosTitle(mySOS)))


###################################################
### code chunk number 44: metadataExtraction1 (eval = FALSE)
###################################################
## sosOfferings(mySOS)
## sosOfferings(mySOS, name = "Rain")


###################################################
### code chunk number 45: metadataExtraction1
###################################################
summary(sosOfferings(mySOS)[[1]])


###################################################
### code chunk number 46: metadataExtraction2 (eval = FALSE)
###################################################
## sosOfferingIds(mySOS)
## names(sosOfferings(mySOS))
## sosName(sosOfferings(mySOS))


###################################################
### code chunk number 47: metadataExtraction3
###################################################
off.temp <- sosOfferings(mySOS)[["ATMOSPHERIC_TEMPERATURE"]]


###################################################
### code chunk number 48: metadataExtraction4
###################################################
off.temp.id <- sosId(off.temp)
off.temp.name <- sosName(off.temp)


###################################################
### code chunk number 49: metadataExtraction4b
###################################################
sosResultModels(off.temp)
sosResponseMode(off.temp)
sosResponseFormats(off.temp)


###################################################
### code chunk number 50: metadataExtraction5
###################################################
off.temp.boundedBy <- sosBoundedBy(off.temp)


###################################################
### code chunk number 51: metadataExtraction6
###################################################
off.temp.boundedBy.bbox <- sosBoundedBy(off.temp, bbox = TRUE)


###################################################
### code chunk number 52: metadataExtraction7
###################################################
off.temp.time <- sosTime(off.temp)
str(off.temp.time)

off.temp.time@beginPosition@time
off.temp.time@endPosition@time
class(off.temp.time@endPosition@time)


###################################################
### code chunk number 53: metadataExtraction8
###################################################
off.temp.time.converted <- sosTime(off.temp, convert = TRUE)
str(off.temp.time.converted)


###################################################
### code chunk number 54: metadataExtraction9
###################################################
sosProcedures(off.temp)
sosObservedProperties(off.temp)
sosFeaturesOfInterest(off.temp)


###################################################
### code chunk number 55: metadataExtraction10
###################################################
sosProcedures(mySOS)[1:2]
sosObservedProperties(mySOS)[1:2]
sosFeaturesOfInterest(mySOS)[1:2]


###################################################
### code chunk number 56: metadataExtraction10
###################################################
sosProcedures(sosOfferings(mySOS)[4:5])
sosObservedProperties(sosOfferings(mySOS)[4:5])
sosFeaturesOfInterest(sosOfferings(mySOS)[3:4])


###################################################
### code chunk number 57: getObservation (eval = FALSE)
###################################################
## getObservation(sos = mySOS, offeringy = myOffering, ...)


###################################################
### code chunk number 58: defaultValue
###################################################
defaultResponseFormatGetObs <- gsub(pattern = "&quot;", replacement = "'", x = sosDefaultGetObsResponseFormat)


###################################################
### code chunk number 59: getObsPropPhen (eval = FALSE)
###################################################
## obs.temp.procedure.1 <- getObservation(sos = mySOS,
## 	offering = off.temp,
## 	procedure = sosProcedures(off.temp)[[2]])
## 
## obs.temp.offering.34 <- getObservation(sos = mySOS,
## 	offering = off.temp,
## 	procedure = sosProcedures(off.temp)[3:4],
## 	observedProperty =
## 		sosObservedProperties(mySOS)[3:4])


###################################################
### code chunk number 60: loadOrDownload_ObsTemp
###################################################
.obsFile <- "obs.temp"
if (.goOnline) {
	obs.temp <- getObservation(sos = mySOS,
		offering = off.temp,
		eventTime = sosCreateTime(sos = mySOS, time = "2009-08-20::2009-08-21"),
		saveOriginal = .obsFile)
} else {
	obs.temp <- parseFile(mySOS, .getFilePath(.obsFile))
}


###################################################
### code chunk number 61: getObs0a (eval = FALSE)
###################################################
## obs.temp <- getObservation(sos = mySOS,
## 		offering = off.temp,
## 		eventTime = sosCreateTime(sos = mySOS, time = "2009-08-20::2009-08-21"),
## 		saveOriginal = .obsFile)


###################################################
### code chunk number 62: getObs0b
###################################################
.results <- lapply(obs.temp, sosResult)
.resultLength <- sapply(.results, nrow)
cat("[sos4R] Received response (size: 27056 bytes), starting parsing ...\n",
	"[sos4R] Finished getObservation to", sosUrl(mySOS), "\n",
	"	--> received", length(obs.temp), "observation(s) having",
	sum(.resultLength), "result values [", toString(.resultLength), "].")


###################################################
### code chunk number 63: getObs1
###################################################
class(obs.temp)
str(obs.temp, max.level = 2)


###################################################
### code chunk number 64: getObs1
###################################################
length(obs.temp)
obs.temp[[1]]
summary(obs.temp)
summary(obs.temp[[1]])


###################################################
### code chunk number 65: getObs1 (eval = FALSE)
###################################################
## obs.temp[2:3]


###################################################
### code chunk number 66: getObs1
###################################################
index.foiId <- sosFeatureIds(obs.temp)[[1]]
index.foiId
obs.temp[index.foiId]

index.obsProp <- sosObservedProperties(off.temp)
obs.temp[index.obsProp]

index.proc <- sosProcedures(obs.temp)[1:4]
index.proc.alternative1 <- sosProcedures(off.temp)[1:4]
index.proc.alternative2 <- sosProcedures(mySOS)
obs.temp[index.proc]


###################################################
### code chunk number 67: getObs2
###################################################
obs.temp.result.2 <- sosResult(obs.temp[[2]])
obs.temp.result <- sosResult(obs.temp[1:2])


###################################################
### code chunk number 68: getObs3
###################################################
temperature.attrs <- attributes(
	obs.temp.result[["urn:ogc:def:property:OGC::Temperature"]])


###################################################
### code chunk number 69: getObsSpatial1
###################################################
obs.temp.foiIDs <- sosFeatureIds(obs.temp)
obs.temp.coords <- sosCoordinates(obs.temp)
obs.temp.coords.1 <- sosCoordinates(obs.temp[[1]])


###################################################
### code chunk number 70: getObsSpatial1
###################################################
sosBoundedBy(obs.temp)
sosBoundedBy(obs.temp, bbox = TRUE)


###################################################
### code chunk number 71: getObsSpatial2
###################################################
result.names <- names(obs.temp.result)
coords.names <- names(obs.temp.coords)
print(toString(result.names))
print(toString(coords.names))

obs.temp.data <- merge(
	x = obs.temp.result,
	y = obs.temp.coords,
	by.x = result.names[[2]],
	by.y = coords.names[[4]])


###################################################
### code chunk number 72: getObsSpatial3
###################################################
obs.temp.data <- merge(x = obs.temp.result,
	y = obs.temp.coords)
str(obs.temp.data, max.level = 2)


###################################################
### code chunk number 73: getObsSpatial4 (eval = FALSE)
###################################################
## sosResult(obs.temp, coordinates = TRUE)


###################################################
### code chunk number 74: temporalFiltering1a
###################################################
# temporal interval creation based on POSIXt classes:
lastWeek.period <- sosCreateTimePeriod(sos = mySOS,
	begin = (Sys.time() - 3600 * 24 * 7),
	end = Sys.time())

oneWeek.period <- sosCreateTimePeriod(sos = mySOS,
		begin = as.POSIXct("2010/01/01"),
		end = as.POSIXct("2010/01/07"))
oneWeek.eventTime <- sosCreateEventTimeList(oneWeek.period)


###################################################
### code chunk number 75: temporalFiltering1b
###################################################
sosCreateTime(sos = mySOS, time = "2007-07-07 07:00::2008-08-08 08:00")
sosCreateTime(sos = mySOS, time = "2007-07-07 07:00/2010-10-10 10:00")

sosCreateTime(sos = mySOS, time = "::2007-08-05")
sosCreateTime(sos = mySOS, time = "2007-08-05/")


###################################################
### code chunk number 76: loadOrDownload_TemporalFiltering2
###################################################
.obsFile <- "obs.oneWeek"
if(.goOnline) {
	obs.oneWeek <- getObservation(sos = mySOS,
		offering = off.temp,
		procedure = sosProcedures(off.temp),
		eventTime = oneWeek.eventTime,
		saveOriginal = .obsFile)
} else {
	obs.oneWeek <- parseFile(mySOS, .getFilePath(.obsFile))
}


###################################################
### code chunk number 77: temporalFiltering2a (eval = FALSE)
###################################################
## obs.oneWeek <- getObservation(sos = mySOS,
## 	offering = off.temp,
## 	# actually not required, as default is 'all procedures':
## 	procedure = sosProcedures(off.temp),
## 	eventTime = oneWeek.eventTime)


###################################################
### code chunk number 78: temporalFiltering2b
###################################################
obs.oneWeek.result <- sosResult(obs.oneWeek)
summary(obs.oneWeek.result[,"urn:ogc:def:property:OGC::Temperature"])


###################################################
### code chunk number 79: temporalFiltering3
###################################################
lastDay.instant <- sosCreateTimeInstant(
	time = as.POSIXct(Sys.time() - 3600 * 24), sos = mySOS)
lastDay.eventTime <- sosCreateEventTime(time = lastDay.instant,
	operator = SosSupportedTemporalOperators()[["TM_After"]])
print(lastDay.eventTime)


###################################################
### code chunk number 80: spatialFiltering1
###################################################
sept09.period <- sosCreateTimePeriod(sos = mySOS,
	begin = as.POSIXct("2009-09-01 00:00"),
	end = as.POSIXct("2009-09-30 00:00"))
sept09.eventTimeList <- sosCreateEventTimeList(
	sept09.period)


###################################################
### code chunk number 81: loadOrDownload_SpatialFiltering
###################################################
.obsFile <- "obs.sept09"
.obsFile2 <- "obs.sept09.bbox"
if(.goOnline) {
	obs.sept09 <- getObservation(sos = mySOS,
		offering = off.temp,
		eventTime = sept09.eventTimeList,
		saveOriginal = .obsFile)

	request.bbox <- sosCreateBBOX(lowLat = 50.0, lowLon = 5.0,
			uppLat = 55.0, uppLon = 10.0,
			srsName = "urn:ogc:def:crs:EPSG:4326")
	request.bbox.foi <- sosCreateFeatureOfInterest(
			spatialOps = request.bbox)
	
	obs.sept09.bbox <- getObservation(sos = mySOS,
			offering = off.temp,
			featureOfInterest = request.bbox.foi,
			eventTime = sept09.eventTimeList,
			saveOriginal = .obsFile2)
} else {
	obs.sept09 <- parseFile(mySOS, .getFilePath(.obsFile))
	obs.sept09.bbox <- parseFile(mySOS, .getFilePath(.obsFile2))
}


###################################################
### code chunk number 82: spatialFiltering2 (eval = FALSE)
###################################################
## obs.sept09 <- getObservation(sos = mySOS,
## 	offering = off.temp,
## 	eventTime = sept09.eventTimeList)


###################################################
### code chunk number 83: spatialFiltering3 (eval = FALSE)
###################################################
## request.bbox <- sosCreateBBOX(lowLat = 50.0, lowLon = 5.0,
## 	uppLat = 55.0, uppLon = 10.0,
## 	srsName = "urn:ogc:def:crs:EPSG:4326")
## request.bbox.foi <- sosCreateFeatureOfInterest(
## 	spatialOps = request.bbox)
## obs.sept09.bbox <- getObservation(sos = mySOS,
## 	offering = off.temp,
## 	featureOfInterest = request.bbox.foi,
## 	eventTime = sept09.eventTimeList)


###################################################
### code chunk number 84: spatialFiltering1
###################################################
print(sosCoordinates(obs.sept09)[,1:2])
print(sosCoordinates(obs.sept09.bbox)[,1:2])


###################################################
### code chunk number 85: featureFiltering1
###################################################
off.temp.fois <- sosFeaturesOfInterest(off.temp)
request.fois <- sosCreateFeatureOfInterest(
	objectIDs = list(off.temp.fois[[1]]))
encodeXML(obj = request.fois, sos = mySOS)


###################################################
### code chunk number 86: loadOrDownload_FeatureFiltering
###################################################
.obsFile <- "obs.oneWeek.fois"
if(.goOnline) {
	obs.oneWeek.fois <- getObservation(sos = mySOS,
		offering = off.temp,
		featureOfInterest = request.fois,
		eventTime = oneWeek.eventTime,
		saveOriginal = .obsFile)
} else {
	obs.oneWeek.fois <- parseFile(mySOS, .getFilePath(.obsFile))
}


###################################################
### code chunk number 87: featureFiltering1a (eval = FALSE)
###################################################
## obs.oneWeek.fois <- getObservation(sos = mySOS,
## 	offering = off.temp,
## 	featureOfInterest = request.fois,
## 	eventTime = oneWeek.eventTime)


###################################################
### code chunk number 88: featureFiltering1b
###################################################
print(sosFeaturesOfInterest(obs.oneWeek.fois))


###################################################
### code chunk number 89: valueFiltering1
###################################################
# result filtering
filter.value <- -2.3
filter.propertyname <- xmlNode(name = ogcPropertyNameName,
	namespace = ogcNamespacePrefix)
xmlValue(filter.propertyname) <-
		"urn:ogc:def:property:OGC::Temperature"
filter.literal <- xmlNode(name = ogcLiteralName,
	namespace = ogcNamespacePrefix)
xmlValue(filter.literal) <- as.character(filter.value)
filter.comparisonop <- xmlNode(
	name = ogcComparisonOpGreaterThanName,
	namespace = ogcNamespacePrefix,
	.children = list(filter.propertyname,
	filter.literal))
filter.result <- xmlNode(name = sosResultName,
	namespace = sosNamespacePrefix,
	.children = list(filter.comparisonop))


###################################################
### code chunk number 90: valueFiltering2
###################################################
filter.result


###################################################
### code chunk number 91: valueFiltering3a (eval = FALSE)
###################################################
## obs.oneWeek <- getObservation(sos = mySOS,
## 	eventTime = oneWeek.eventTime,
## 	offering = sosOfferings(mySOS)[["ATMOSPHERIC_TEMPERATURE"]])


###################################################
### code chunk number 92: loadOrDownload_ValueFiltering
###################################################
.obsFile <- "obs.oneWeek.filter"
if(.goOnline) {
	obs.oneWeek.filter <- getObservation(sos = mySOS,
		eventTime = oneWeek.eventTime,
		offering = sosOfferings(mySOS)[["ATMOSPHERIC_TEMPERATURE"]],
		result = filter.result,
		saveOriginal = .obsFile)
} else {
	obs.oneWeek.filter <- parseFile(mySOS, .getFilePath(.obsFile))
}


###################################################
### code chunk number 93: valueFiltering3b (eval = FALSE)
###################################################
## # request  values for the week with a value higher than 0 degrees:
## obs.oneWeek.filter <- getObservation(sos = mySOS,
## 	eventTime = oneWeek.eventTime,
## 	offering = sosOfferings(mySOS)[["ATMOSPHERIC_TEMPERATURE"]],
## 	result = filter.result)


###################################################
### code chunk number 94: valueFiltering3c
###################################################
print(paste("Filtered:", dim(sosResult(obs.oneWeek.filter))[[1]], 
	"-vs.- Unfiltered:", dim(sosResult(obs.oneWeek))[[1]]))


###################################################
### code chunk number 95: resultExporting1
###################################################
library("sp")


###################################################
### code chunk number 96: resultExporting2 (eval = FALSE)
###################################################
## obs.oneWeek <- getObservation(sos = mySOS,
## 	offering = off.temp,
## 	procedure = sosProcedures(off.temp),
## 	eventTime = oneWeek.eventTime)


###################################################
### code chunk number 97: resultExporting3
###################################################
# Create SpatialPointsDataFrame from result features
coords <- sosCoordinates(obs.oneWeek[[1]])
crs <- sosGetCRS(obs.oneWeek[[1]])
spdf <- SpatialPointsDataFrame(coords = coords[,1:2],
	data = data.frame(coords[,4]), proj4string = crs)
str(spdf)


###################################################
### code chunk number 98: getObsCRS
###################################################
sosGetCRS(obs.temp)
sosGetCRS(obs.oneWeek)


###################################################
### code chunk number 99: loadOrDownload_GetObsById
###################################################
.obsFile <- "obsId"
if(.goOnline) {
	obsId <- getObservationById(sos = mySOS,
		observationId = "o_3508493",
		saveOriginal = .obsFile)
} else {
	obsId <- parseFile(mySOS, .getFilePath(.obsFile))
}


###################################################
### code chunk number 100: getObsById1a (eval = FALSE)
###################################################
## obsId <- getObservationById(sos = mySOS,
## 	observationId = "o_3508493")


###################################################
### code chunk number 101: getObsById1b
###################################################
sosResult(obsId, coordinates = TRUE)


###################################################
### code chunk number 102: getObsById3 (eval = FALSE)
###################################################
## # generated file name, find file in working directory:
## obsId <- getObservationById(sos = mySOS,
## 	observationId = "o_3508493",
## 	saveOriginal = TRUE)
## .files <- list.files(getwd())
## .observationFiles <- c()
## for(.f in .files) { # %in% not working with Sweave
## 	if(length(grep("^o_", .f, value=TRUE)) > 0)
## 		.observationFiles <- c(.observationFiles, .f)
## }
## obsId <- parseFile(sos = mySOS,
## 	file = .observationFiles[[1]])
## 


###################################################
### code chunk number 103: getObsById3 (eval = FALSE)
###################################################
## # manually selected file name:
## obsId <- getObservationById(sos = mySOS,
## 	#verbose = TRUE,
## 	observationId = "o_3508493",
## 	saveOriginal = "myObservation")


###################################################
### code chunk number 104: inclusionExclusion
###################################################
parsers <- SosParsingFunctions(
	"ExceptionReport" = function() {
		return("Got Exception!")
	},
	include = c("GetObservation", "ExceptionReport"))
print(names(parsers))

parsers <- SosParsingFunctions(
		"ExceptionReport" = function() {
			return("Got Exception!")
		},
		include = c("GetCapabilities"))
print(names(parsers))


###################################################
### code chunk number 105: encoders1 (eval = FALSE)
###################################################
## sosEncoders(mySOS)


###################################################
### code chunk number 106: encoders2
###################################################
names(sosEncoders(mySOS))


###################################################
### code chunk number 107: encoders3 (eval = FALSE)
###################################################
## myPostEncoding <- function(object, sos, verbose) {
## 	return(str(object))
## }
## # Connection will not be establihsed because of mising objects
## mySOS2 = SOS(sosUrl(mySOS),
## 	encoders = SosEncodingFunctions("POST" = myPostEncoding))


###################################################
### code chunk number 108: encoders4
###################################################
showMethods("encodeXML")
showMethods("encodeKVP")


###################################################
### code chunk number 109: encoders5 (eval = FALSE)
###################################################
## setMethod(f = "encodeXML",
##   signature = signature(obj = "POSIXt", sos = "SOS"),
##     def = function(obj, sos, verbose) {
##       if(verbose) cat("Using my own time encoding... ")
## 
##       # time zone hack to fix that the time format option
##       # %z does not work on windows machines:
##       .time <- obj + 11 * 60 * 60 # add 11 hours
##       .formatted <- strftime(x = .time,
##         format = sosTimeFormat(sos))
##       .formatted <- paste(.formatted, 
##         "+11:00", sep = "")	# append 11:00
## 
##       if(verbose) cat("Formatted ", toString(obj),
##         " to ", .formatted, "\n")
##       return(.formatted)
##     }
## )


###################################################
### code chunk number 110: parsers1 (eval = FALSE)
###################################################
## sosParsers(mySOS)


###################################################
### code chunk number 111: parsers2
###################################################
names(sosParsers(mySOS))


###################################################
### code chunk number 112: loadOrDownload_Parsers3
###################################################
.obsFile <- "err.response.RData"
.sos2 <- "mySOS2.RData"
if(.goOnline) {
	myER <- function(xml) {
		return("EXCEPTION!!!11")
	}
	myParsers <- SosParsingFunctions("ExceptionReport" = myER)
	mySOS2 <- SOS(sosUrl(mySOS), parsers = myParsers)
	save(mySOS2, file = .sos2)
	
	err.response <- getObservation(mySOS2,
		offering = sosOfferings(mySOS2)[[1]],
		observedProperty = list("Bazinga"))
	save(err.response, file = .obsFile)
} else {
	load(.sos2)
	load(.obsFile)
}


###################################################
### code chunk number 113: parsers3 (eval = FALSE)
###################################################
## # Create own parsing function:
## myER <- function(xml) {
## 	return("EXCEPTION!!!11")
## }
## myParsers <- SosParsingFunctions("ExceptionReport" = myER)
## mySOS2 <- SOS(sosUrl(mySOS), parsers = myParsers)
## # Triggers exception:
## err.response <- getObservation(mySOS2, verbose = TRUE,
## 	offering = sosOfferings(mySOS2)[[1]],
## 	observedProperty = list("Bazinga!"))


###################################################
### code chunk number 114: parsers3
###################################################
print(err.response)


###################################################
### code chunk number 115: parsers4a (eval = FALSE)
###################################################
## SosDisabledParsers()


###################################################
### code chunk number 116: parsers4b
###################################################
names(SosDisabledParsers())


###################################################
### code chunk number 117: loadOrDownload_Parsers5
###################################################
.obsFile <- "response.noparsing.RData"
.sos2.dis <- "mySOS2.disabled.RData"
if(.goOnline) {
	mySOS2.disabled <- SOS(sosUrl(mySOS), # verboseOutput = TRUE,
			parsers = SosDisabledParsers())
	response.noparsing <- getObservation(mySOS2.disabled,
			offering = sosOfferings(mySOS2.disabled)[[1]],
			observedProperty = list("Bazinga!"))
	save(mySOS2.disabled, file = .sos2.dis)
#	save(response.noparsing, file = .obsFile)
} else {
	load(.sos2.dis)
#	load(.obsFile)
}


###################################################
### code chunk number 118: parsers5a (eval = FALSE)
###################################################
## mySOS2.disabled <- SOS(sosUrl(mySOS),
## 		parsers = SosDisabledParsers())
## response.noparsing <- getObservation(mySOS2.disabled,
## 	offering = sosOfferings(mySOS2.disabled)[[1]],
## 	observedProperty = list("Bazinga"))


###################################################
### code chunk number 119: parsers5b (eval = FALSE)
###################################################
## class(response.noparsing)
## print(xmlName(xmlRoot(response.noparsing)))
## # (Using XML functions here for accesing the root of a
## # document and the name of an element.)


###################################################
### code chunk number 120: converters0
###################################################
value <- 2.0
value.string <- sosConvertString(x = value, sos = mySOS)
print(class(value.string))

value <- "2.0"
value.double <- sosConvertDouble(x = value, sos = mySOS)
print(class(value.double))

value <- "1"
value.logical <- sosConvertLogical(x = value, sos = mySOS)
print(class(value.logical))

value <- "2010-01-01T12:00:00.000"
value.time <- sosConvertTime(x = value, sos = mySOS)
print(class(value.time))


###################################################
### code chunk number 121: converters1
###################################################
names(SosDataFieldConvertingFunctions())


###################################################
### code chunk number 122: converters2 (eval = FALSE)
###################################################
## sosDataFieldConverters(mySOS)


###################################################
### code chunk number 123: loadOrDownload_Converters3
###################################################
.obsFile <- "mbariObs1"
.mbari <- "mbari.RData"
if(.goOnline) {
	MBARI <- SOS("http://mmisw.org/oostethys/sos",
		method = SosSupportedConnectionMethods()[["GET"]])
	myOff <- sosOfferings(MBARI)[[1]]
	myProc <- sosProcedures(MBARI)[[1]]

	mbariObs1 <- try(
			getObservation(sos = MBARI, offering = myOff, verbose = TRUE,
				procedure = myProc, responseFormat = NA_character_,
				saveOriginal = .obsFile)
		)
	save(MBARI, file = .mbari)
} else {
	load(.mbari)
	mbariObs1 <- parseFile(MBARI, .getFilePath(.obsFile))
}


###################################################
### code chunk number 124: converters3a (eval = FALSE)
###################################################
## # GET connection
## MBARI <- SOS("http://mmisw.org/oostethys/sos",
## 	method = SosSupportedConnectionMethods()[["GET"]])
## myOff <- sosOfferings(MBARI)[[1]]
## myProc <- sosProcedures(MBARI)[[1]]
## mbariObs1 <- try(
## 	getObservation(sos = MBARI, offering = myOff,
## 		procedure = myProc, responseFormat = NA_character_)
## )


###################################################
### code chunk number 125: converters3b (eval = FALSE)
###################################################
## warnings()


###################################################
### code chunk number 126: converters3c
###################################################
# ran manually:
#converterWarnings <- warnings()[25:29]
#.convWarnFile <- "converterWarnings.Rdata"
#save(converterWarnings, file = .convWarnFile)
#load(.convWarnFile)
#print(converterWarnings)
cat("...")
cat("25: In FUN(X[[7L]], ...) :
		swe:Quantity given without unit of measurement: Salinity
26: In .valParser(values = obj[[sweValuesName]], fields = .fields,  ... :
				No converter for the unit of measurement  S/m  with the definition  http://mmisw.org/ont/cf/parameter/conductivity ! Trying a default, but you can add one when creating a SOS using SosDataFieldConvertingFunctions().
27: In .valParser(values = obj[[sweValuesName]], fields = .fields,  ... :
				No converter found! Skipping field Conductivity 
No converter found! Skipping field http://mmisw.org/ont/cf/parameter/conductivity 
No converter found! Skipping field S/m 

28: In .valParser(values = obj[[sweValuesName]], fields = .fields,  ... :
				No converter found for the given field Salinity, http://mmisw.org/ont/cf/parameter/sea_water_salinity, NA
29: In .valParser(values = obj[[sweValuesName]], fields = .fields,  ... :
				No converter found! Skipping field Salinity 
No converter found! Skipping field http://mmisw.org/ont/cf/parameter/sea_water_salinity 
No converter found! Skipping field NA 

30: In FUN(X[[7L]], ...) :
		swe:Quantity given without unit of measurement: Salinity")
cat("...")


###################################################
### code chunk number 127: converters4
###################################################
myConverters <- SosDataFieldConvertingFunctions(
	"S/m" = sosConvertDouble,
	"http://mmisw.org/ont/cf/parameter/sea_water_salinity"
			= sosConvertDouble)


###################################################
### code chunk number 128: loadOrDownload_Converters4
###################################################
.obsFile <- "mbariObs2"
.mbari2 <- "mbari2.RData"
if(.goOnline) {
	MBARI2 <- SOS("http://mmisw.org/oostethys/sos",
		method = SosSupportedConnectionMethods()[["GET"]],
		dataFieldConverters = myConverters)
	mbariObs2 <- getObservation(sos = MBARI2, offering = myOff, #verbose = TRUE,
		procedure = myProc, responseFormat = NA_character_,
		saveOriginal = .obsFile)
	save(MBARI2, file = .mbari2)
} else {
	load(.mbari2)
	mbariObs2 <- parseFile(MBARI2, .getFilePath(.obsFile))
}


###################################################
### code chunk number 129: converters4a (eval = FALSE)
###################################################
## MBARI2 <- SOS("http://mmisw.org/oostethys/sos",
## 	method = SosSupportedConnectionMethods()[["GET"]],
## 	dataFieldConverters = myConverters)
## mbariObs2 <- getObservation(sos = MBARI2, offering = myOff,
## 	procedure = myProc, responseFormat = NA_character_)


###################################################
### code chunk number 130: converters5
###################################################
toString(names(sosResult(mbariObs1)))
toString(names(sosResult(mbariObs2)))


###################################################
### code chunk number 131: exceptionData (eval = FALSE)
###################################################
## OwsExceptionsData()


###################################################
### code chunk number 132: exceptionTable
###################################################
library(xtable)
print(xtable(x = OwsExceptionsData()[,1:3],
				caption = "Exception Data Table (without HTTP columns).",
				label = c("tab:exceptions"), table.placement = "tbp",
				align = c("l", "l", "p{5cm}", "p{3cm}"),
				caption.placement = "top"),
		include.rownames = FALSE)


###################################################
### code chunk number 133: exceptionWarning1 (eval = FALSE)
###################################################
## response <- try(getObservationById(sos = mySOS,
## 	observationId = "o_not_there"))


###################################################
### code chunk number 134: exceptionWarning2b
###################################################
cat("Warning:
In .handleExceptionReport(sos, .response) :
  Object of class OwsExceptionReport; version: 1.0.0; lang: NA;
 1 exception(s) (code @ locator : text):
  NoApplicableCode @ NA :
	Error while creating observations from database query result set: ERROR: invalid input syntax for integer: \"not_there\" ")


###################################################
### code chunk number 135: exceptionWarning2c
###################################################
.file <- "exceptionResponse.RData"
if(.goOnline) {
	response <- getObservationById(sos = mySOS,
					observationId = "o_not_there")
	save(response, file = .file)
} else {
	load(.file)
}


###################################################
### code chunk number 136: exceptionWarning2d
###################################################
response


###################################################
### code chunk number 137: inspect (eval = FALSE)
###################################################
## off.4 <- sosOfferings(mySOS)[[4]]
## getObservation(sos = mySOS, offering = off.4,
## 	procedure = sosProcedures(off.4)[[1]],
## 	inspect = TRUE)
## 
## getObservation(sos = mySOS,	offering = off.4,
## 	procedure = sosProcedures(off.4)[[1]],
## 	vebose = TRUE)
## 
## verboseSOS <- SOS(sosUrl(mySOS), verboseOutput = TRUE)


###################################################
### code chunk number 138: install (eval = FALSE)
###################################################
## install.packages("sos4R")


###################################################
### code chunk number 139: demo (eval = FALSE)
###################################################
## # list available demos:
## demo(package = "sos4R")
## 
## # run a demo:
## demo("airquality")


###################################################
### code chunk number 140: exampleServices
###################################################
SosExampleServices()


###################################################
### code chunk number 141: cheatSheet (eval = FALSE)
###################################################
## sosCheatSheet()


###################################################
### code chunk number 142: options
###################################################
print("###################################################################")
if(.goOnline) {
	print("####### Data for this vignette was downloaded ONLINE. #############")
} else {
	print("####### Data for this vignette was loaded from LOCAL FILES. #######")
}
print("###################################################################")


