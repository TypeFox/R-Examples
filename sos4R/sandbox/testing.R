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
# the terms of the GNU General Publipc License version 2 as published by the   #
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
# Created: 2010-06-18                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

################################################################################
# This file is the main development basis for the package sos4R. Examples should
# eventually be moved to a Vignette.
################################################################################

################################################################################
source("/home/daniel/Dropbox/2010_SOS4R/workspace/sos4R/sandbox/loadSources.R")


################################################################################
# RCurl

sosUrl = "http://giv-sos.uni-muenster.de:8080/ClimateSOS/sos"
request = "service=SOS&request=GetCapabilities&acceptVersions=1.0.0,2.0.0&sections=OperationsMetadata,ServiceIdentification,ServiceProvider,Filter_Capabilities,Contents&acceptFormats=text/xml"
url = paste(sosUrl, request, sep = "?")
getURL(url, verbose = TRUE)


################################################################################
# GetCapabilities

climatesosUrl = "http://giv-sos.uni-muenster.de:8080/ClimateSOS/sos"
climatesos = SOS(climatesosUrl)
weathersosUrl = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos"
weathersos = SOS(weathersosUrl, verboseOutput = FALSE)

caps = getCapabilities(weathersos)
sosCaps(weathersos)


################################################################################
# DescribeSensor

climatesosUrl = "http://giv-sos.uni-muenster.de:8080/ClimateSOS/sos"
climatesos = SOS(climatesosUrl, verboseOutput = FALSE)
id = "urn:ogc:object:feature:WMOStation:10280"
describeSensor(sos = climatesos, procedure = id, saveOriginal = "D:/text.xml")

# !!! describeSensor does not check if using GET, because Capabilities lack that DCP in current SOS!
weathersosUrl = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos"
weathersos = SOS(weathersosUrl, method = "POST", verboseOutput = TRUE)
sensor <- describeSensor(weathersos, sosProcedures(weathersos)[[1]])
sensor <- describeSensor(weathersos, "manniK")
sensor <- describeSensor(weathersos, sosProcedures(weathersos)[[1]], verbose = FALSE)

################################################################################
# GetObservation

# mandatory:
go.service = "SOS"
go.version = "1.0.0"
go.offering = "region_3"
go.observedProperty = c("urn:ogc:def:property:OGC:1.0:temperature", "urn:ogc:def:property:OGC:1.0:windDirection")
go.responseFormat = "text/xml;subtype=&quot;om/1.0.0&quot;"

# creation method
go <- SosGetObservation(service = go.service, version = go.version, offering = go.offering, observedProperty =  go.observedProperty, responseFormat =  go.responseFormat)
go.xml <- encode(go)

################################################################################
# ExceptionReports

weathersos.url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos"
weathersos = SOS(weathersos.url)

id.correct = "urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93"
id.incorrect = "lala"

er.xmltext <- "<ows:ExceptionReport xmlns:ows=\"http://www.opengis.net/ows/1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"1.0.0\" xsi:schemaLocation=\"http://schemas.opengis.net/ows/1.1.0/owsExceptionReport.xsd\"><ows:Exception exceptionCode=\"VersionNegotiationFailed\" locator=\"AcceptVersions\"><ows:ExceptionText>The parameter 'AcceptVersions' does not contain the version of this SOS: '1.0.0'</ows:ExceptionText></ows:Exception></ows:ExceptionReport>"
er.xmltext2 <- "<ows:ExceptionReport xmlns:ows=\"http://www.opengis.net/ows/1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"1.0.0\" lang=\"de-DE\" xsi:schemaLocation=\"http://schemas.opengis.net/ows/1.1.0/owsExceptionReport.xsd\"><ows:Exception exceptionCode=\"VersionNegotiationFailed\" locator=\"AcceptVersions\"><ows:ExceptionText>The parameter 'AcceptVersions' does not contain the version of this SOS: '1.0.0'</ows:ExceptionText></ows:Exception><ows:Exception exceptionCode=\"NoApplicableCode\" locator=\"@home\"><ows:ExceptionText>Just a second exception to make things saver...</ows:ExceptionText></ows:Exception></ows:ExceptionReport>"
er.doc <- xmlParseDoc(er.xmltext)
er.root <- xmlRoot(er.doc)

er2.doc <- xmlParseDoc(er.xmltext2)
er2.parsed <- parseOwsExceptionReport(er2.doc)
str(er2.parsed)

# unparsable report for wrong getcap:
er.xmltext3 <- "<ows:ExceptionReport xmlns:ows=\"http://www.opengis.net/ows/1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"1.0.0\" xsi:schemaLocation=\"http://schemas.opengis.net/ows/1.1.0/owsExceptionReport.xsd\"><ows:Exception exceptionCode=\"VersionNegotiationFailed\" locator=\"AcceptVersions\"><ows:ExceptionText>The parameter 'AcceptVersions' HAR HAR does not contain the version of this SOS: '1.0.0'</ows:ExceptionText></ows:Exception></ows:ExceptionReport>"
str(parseOwsExceptionReport(xmgetlParseDoc(er.xmltext3)))
# works from here!

debug(parseOwsExceptionReport)
debug(parseOwsException)

# no exception text
er.xmltext4 <- "<ows:ExceptionReport xmlns:ows=\"http://www.opengis.net/ows/1.1\" lang=\"de-DE\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"1.0.0\" lang=\"de-DE\" xsi:schemaLocation=\"http://schemas.opengis.net/ows/1.1.0/owsExceptionReport.xsd\" />"


################################################################################
# escaping special characters in url
input1 <- "urn:ogc:def:property:OGC:1.0:temperature_#_1"
input2 <- "om:SpatialObservation"
input3 <- "text/xml;subtype=\"OM/1.0.0\""
input4 <- "urn@lala@home+home+home+lat,lon,lat,lon"

.kvpEscapeSpecialCharacters(input1)
.kvpEscapeSpecialCharacters(input2)
.kvpEscapeSpecialCharacters(input3)
.kvpEscapeSpecialCharacters(input4)

weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos",
		method = SosSupportedConnectionMethods()[[1]], verboseOutput = TRUE)

describeSensor(sos = weathersos, procedure = "sos@home+street:name,one,two")
getObservation(sos = weathersos, offering = sosOfferingIds(weathersos)[[1]])

################################################################################
source("/home/daniel/Dropbox/2010_SOS4R/workspace/sos4R/sandbox/loadSources.R")
################################################################################

################################################################################
# problem with "/" character
url = "http://localhost:8080/ClimateSOS-local/sos"
request1 = "service=SOS&request=GetCapabilities&acceptVersions=1.0.0,2.0.0&sections=All&acceptFormats=text/xml"
request2 = "service=SOS&request=GetCapabilities&acceptVersions=1.0.0,2.0.0&sections=All&acceptFormats=text%2Fxml"

getURL(paste(url, request1, sep = "?"))
getURL(paste(url, request2, sep = "?"))


################################################################################
# Parsing the capabilities file...

weathersos = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
#caps <- getCapabilities(weathersos, verbose = TRUE)
weathersos@capabilities

# some error
caps <- '...'
xmlCaps <- xmlParseDoc(caps)
debug(parseSosCapabilities)
parsedCaps <- parseSosCapabilities(xmlCaps, weathersos)

opXml <- xmlParseDoc('COPY HERE')
debug(parseOwsOperation)
parseOwsOperation(xmlRoot(opXml))


################################################################################
# accessor functions
weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
sosUrl(weathersos)
sosMethod(weathersos)
sosVersion(weathersos)
sosCaps(weathersos)

################################################################################
# POST

# manually:
getCapRequest <- '<?xml version="1.0" encoding="UTF-8"?><GetCapabilities xmlns="http://www.opengis.net/sos/1.0" xmlns:ows="http://www.opengis.net/ows/1.1" xmlns:ogc="http://www.opengis.net/ogc" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xsi:schemaLocation="http://www.opengis.net/sos/1.0 http://schemas.opengis.net/sos/1.0.0/sosGetCapabilities.xsd" service="SOS"><ows:AcceptVersions><ows:Version>1.0.0</ows:Version></ows:AcceptVersions><ows:Sections><ows:Section>OperationsMetadata</ows:Section><ows:Section>ServiceIdentification</ows:Section><ows:Section>ServiceProvider</ows:Section><ows:Section>Filter_Capabilities</ows:Section><ows:Section>Contents</ows:Section></ows:Sections></GetCapabilities>'
# using 'post' for application/x-www-form-urlencoded content
caps.response <- postForm(uri = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos",
		request = getCapRequest,
		style = "POST",
		.encoding = "UTF-8")
caps.doc <- xmlParseDoc(caps.response)
caps <- parseSosCapabilities(caps.doc)

# GetCapabilities
sos = SOS("http://localhost:8080/ClimateSOS-local/sos", method = "POST",
		verboseOutput = TRUE)
sos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos", method = "POST",
		verboseOutput = TRUE)
caps = sosCaps(sos)

# DescribeSensor
procedures = sosProcedures(sos)
sensor.10 <- describeSensor(sos = sos, procedure = procedures[[1]])

################################################################################
# testing to call functions from a list
myFunc1 <- function(xml) {
	print("myfunc: ")
	print(xml)
	return(list(xml, "anotherone"))
}
value <- "lala la"
temp <- list(func1 = myFunc1)
# call the function:
result <- temp[["func1"]](value)
result

################################################################################
# Replace a parsing function... and inspect the request and response
myParseSensorML <- function(obj) {
	root <- xmlRoot(obj)
	return(xmlName(root))
}
myER <- function(xml) {
	return("EXCEPTION! RUN!!!!1111")
}
# testing:
myParsers <- SosParsingFunctions("DescribeSensor" = myParseSensorML)
SosParsingFunctions("DescribeSensor" = myParseSensorML, "ExceptionReport" = myER)
SosParsingFunctions("DescribeSensor" = myParseSensorML, include = c("GetObservation", "DescribeSensor"))
SosParsingFunctions("DescribeSensor" = myParseSensorML, exclude = c("GetObservation", "DescribeSensor"))


weathersos = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos",
		method = "POST",
		parsers = SosParsingFunctions("DescribeSensor" = myParseSensorML, "ExceptionReport" = myER),
		verboseOutput = FALSE)
sensor <- describeSensor(weathersos, sosProcedures(weathersos)[[1]], inspect = TRUE)
# WORKS! YEAH!

################################################################################
# Replace an encoding function
myPostEncoding <- function(object, v) {
	# myPostEncoding
	.request <- encodeRequestXML(obj = object, verbose = v)
	# attach comment node to show this is actually used
	.request[[xmlSize(.request)+1]] <- xmlCommentNode("hej hej!")
	return(.request)
}

weathersos = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos",
		method = "POST",
		encoders = SosEncodingFunctions("POST" = myPostEncoding),
		verboseOutput = TRUE)
sensor <- describeSensor(weathersos, sosProcedures(weathersos)[[1]],
		inspect = TRUE)
# works!

################################################################################
# inspecting XML using dummy parsing function
weathersos2 = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos",
		method = "POST",
		parsers = SosDisabledParsers,
		verboseOutput = FALSE)
sensor2 <- describeSensor(weathersos2, sosProcedures(weathersos2)[[1]])
sensor2 <- describeSensor(weathersos2, "lala")
class(sensor) # from example above
class(sensor2)
# works!

################################################################################
# parsing om:Measurement

meas <- '<?xml version="1.0" encoding="UTF-8"?><om:ObservationCollection xmlns:om="http://www.opengis.net/om/1.0" xmlns:gml="http://www.opengis.net/gml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:sa="http://www.opengis.net/sampling/1.0" gml:id="oc_0" xsi:schemaLocation="http://www.opengis.net/sos/1.0 http://schemas.opengis.net/sos/1.0.0/sosAll.xsd"><om:member><om:Measurement gml:id="o_3376580"><om:samplingTime><gml:TimeInstant xsi:type="gml:TimeInstantType"><gml:timePosition>2010-09-08T09:45:00.000+02:00</gml:timePosition></gml:TimeInstant></om:samplingTime><om:procedure xlink:href="urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93"/><om:observedProperty xlink:href="urn:ogc:def:property:OGC::Temperature"/><om:featureOfInterest><sa:SamplingPoint gml:id="urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93"><gml:name>NOT_SET</gml:name><sa:position><gml:Point><gml:pos srsName="urn:ogc:def:crs:EPSG:4326">51.9412 7.6103</gml:pos></gml:Point></sa:position></sa:SamplingPoint></om:featureOfInterest><om:result uom="Cel">12.9</om:result></om:Measurement></om:member><om:member><om:Measurement gml:id="o_3376580"><om:samplingTime><gml:TimeInstant xsi:type="gml:TimeInstantType"><gml:timePosition>2010-09-08T09:45:00.000+02:00</gml:timePosition></gml:TimeInstant></om:samplingTime><om:procedure xlink:href="urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93"/><om:observedProperty xlink:href="urn:ogc:def:property:OGC::Temperature"/><om:featureOfInterest><sa:SamplingPoint gml:id="urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93"><gml:name>NOT_SET</gml:name><sa:position><gml:Point><gml:pos srsName="urn:ogc:def:crs:EPSG:4326">51.9412 7.6103</gml:pos></gml:Point></sa:position></sa:SamplingPoint></om:featureOfInterest><om:result uom="Cel">12.9</om:result></om:Measurement></om:member></om:ObservationCollection>'
measDoc <- xmlParseDoc(meas)

# for partial testing
tempM <- xmlChildren(xmlRoot(measDoc))[["member"]][["Measurement"]]
parsedM <- parseMeasurement(tempM, SosParsingFunctions(), verbose = T)
tempSP <- xmlChildren(tempM[["featureOfInterest"]])[[1]]
parseSamplingPoint(tempSP)

om <- parseOM(obj = measDoc, parsers = SosParsingFunctions(), verbose = TRUE)
str(om)
# works!

################################################################################
# test handling two of a few not supported observation specializations:
# category observation and geometry observation
obsDoc2 <- xmlParseDoc('<?xml version="1.0" encoding="UTF-8"?><om:ObservationCollection xmlns:om="http://www.opengis.net/om/1.0" xmlns:gml="http://www.opengis.net/gml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:swe="http://www.opengis.net/swe/1.0.1" gml:id="oc_0" xsi:schemaLocation="http://www.opengis.net/sos/1.0 http://schemas.opengis.net/sos/1.0.0/sosAll.xsd"><om:member><om:CategoryObservation></om:CategoryObservation></om:member><om:member><om:GeometryObservation></om:GeometryObservation></om:member></om:ObservationCollection>')
tempObs2 <- parseOM(obsDoc2, parsers = SosParsingFunctions(), verbose = TRUE)
# works!


################################################################################
# creating data frames with time classes
value1 <- list(as.POSIXct("2008-01-01"), as.POSIXct("2009-02-02"),
		as.POSIXct("2010-03-03"))
value2 <- c("lala", "nana", "pooh")
value3 <- c(10.1, 12.4, 17.42)
value4 <- c("10.1", "12.4", "17.42")
value5 <- c(strptime("2010-03-01T12:15:00.000+01:00", sosDefaultTimeFormat),
		strptime("2010-03-02T12:30:00.000+01:00", sosDefaultTimeFormat))
value6 <- c("2010-03-01T12:15:00.000+01:00", "2010-03-02T12:30:00.000+01:00")

df <- data.frame(value1, value2, value3)
df1 <- data.frame(time = value5, name = value2[2:3])
df2 <- data.frame(temperature = value3)
df2[,"lala"] <- values4

timeDF <- data.frame(value5, stringsAsFactors = FALSE)
str(timeDF)
# WORKS, but only if value5 is NOT A LIST

# procedure as in parsing method:
weathersos = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos", method = "POST", verboseOutput = FALSE)
rawTestCurrentValues <- list("2010-03-01T12:15:00.000+01:00",
		"2010-03-02T12:30:00.000+01:00")
testCurrentValues <- lapply(rawTestCurrentValues, sosConvertTime, sos = weathersos)
str(testCurrentValues)
str(data.frame(I(testCurrentValues)))

# looks good... but columnList <- list(testCurrentValues, value2[2:3]) # does not do the trick
cbind(tempTimeDF, tempValuesDF)
str(cbind(tempTimeDF, tempValuesDF))
# why is this not like the following?
str(df1)

str(sosConvertTime(x = rawTestCurrentValues, sos = weathersos))
str(as.double(value4))
# not bad, so just give whole list to conversion function instead of calling lapply
str(sosConvertTime(rawTestCurrentValues, sos = weathersos))

tempTimeDF <- data.frame(sosConvertTime(
				x = list("2010-03-01T12:15:00.000+01:00",
						"2010-03-02T12:30:00.000+01:00"),
				sos = weathersos))
# cannot set name with variable on creation of the data.frame
names(tempTimeDF) <- "time"
tempValuesDF <- data.frame("values" = sosConvertDouble(list("11.1", "12.4")))
tempId = "tempId"
tempData = data.frame(tempId = seq(1:2))
tempData <- cbind(tempId, tempTimeDF, tempValuesDF); str(tempData)
tempData[,!colnames(tempData)%in%tempId]
# YEAH! WORKS!


################################################################################
# parsing om:Observation

obsDoc <- xmlParseDoc('<?xml version="1.0" encoding="UTF-8"?><om:ObservationCollection xmlns:om="http://www.opengis.net/om/1.0" xmlns:gml="http://www.opengis.net/gml" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:xlink="http://www.w3.org/1999/xlink" xmlns:swe="http://www.opengis.net/swe/1.0.1" gml:id="oc_0" xsi:schemaLocation="http://www.opengis.net/sos/1.0 http://schemas.opengis.net/sos/1.0.0/sosAll.xsd"><om:member><om:Observation><om:samplingTime><gml:TimePeriod xsi:type="gml:TimePeriodType"><gml:beginPosition>2010-09-08T09:45:00.000+02:00</gml:beginPosition><gml:endPosition>2010-09-08T09:45:00.000+02:00</gml:endPosition></gml:TimePeriod></om:samplingTime><om:procedure xlink:href="urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93"/><om:observedProperty><swe:CompositePhenomenon gml:id="cpid0" dimension="1"><gml:name>resultComponents</gml:name><swe:component xlink:href="urn:ogc:data:time:iso8601"/><swe:component xlink:href="urn:ogc:def:property:OGC::Temperature"/></swe:CompositePhenomenon></om:observedProperty><om:featureOfInterest><gml:FeatureCollection/></om:featureOfInterest><om:result><swe:DataArray><swe:elementCount><swe:Count><swe:value>1</swe:value></swe:Count></swe:elementCount><swe:elementType name="Components"><swe:SimpleDataRecord><swe:field name="Time"><swe:Time definition="urn:ogc:data:time:iso8601"/></swe:field><swe:field name="feature"><swe:Text definition="urn:ogc:data:feature"/></swe:field><swe:field name="urn:ogc:def:property:OGC::Temperature"><swe:Quantity definition="urn:ogc:def:property:OGC::Temperature"><swe:uom code="Cel"/></swe:Quantity></swe:field></swe:SimpleDataRecord></swe:elementType><swe:encoding><swe:TextBlock decimalSeparator="." tokenSeparator="," blockSeparator=";"/></swe:encoding><swe:values>2010-03-01T12:15:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,73.0;2010-03-01T12:30:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,68.0;2010-03-01T12:45:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,69.0;2010-03-01T13:00:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,65.0;2010-03-01T13:15:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,61.0;2010-03-01T13:30:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,56.0;2010-03-01T13:45:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,60.0;2010-03-01T14:00:00.000+01:00,urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93,57.0;</swe:values></swe:DataArray></om:result></om:Observation></om:member></om:ObservationCollection>')
weathersos = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos", method = "POST", verboseOutput = FALSE)

# testing parts of the observation document:
tempPP <- parsePhenomenonProperty(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omObservedPropertyName]])
tempCP <- parseCompositePhenomenon(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omObservedPropertyName]][[sweCompositePhenomenonName]])
tempTP <- parseSamplingTime(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omSamplingTimeName]])
tempFOI <- parseFOI(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omFeatureOfInterestName]])
tempResult <- parseResult(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omResultName]], sos = weathersos, verbose = TRUE)
tempDA <- parseDataArray(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omResultName]][[sweDataArrayName]], sos = weathersos, verbose = TRUE)
tempFields <- parseElementType(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omResultName]][[sweDataArrayName]][[sweElementTypeName]])
tempEncoding <- parseEncoding(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omResultName]][[sweDataArrayName]][[sweEncodingName]])
tempValues <- parseValues(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]][[omResultName]][[sweDataArrayName]][[sweValuesName]], fields = tempFields, encoding = tempEncoding, sos = weathersos, verbose = TRUE)

tempObs <- parseObservation(xmlRoot(obsDoc)[[omMemberName]][[omObservationName]], sos = weathersos, verbose = TRUE)

# if the observation has a foi with a coordinate, add that to all data rows?

################################################################################
source("/home/daniel/Dokumente/2010_SOS4R/workspace/sos4R/sandbox/loadSources.R")
################################################################################
# GetObservationById

# some ids for weathersos
ids = paste("o_33765", c(80:99), sep = "")
getobsbyid <- GetObservationById("SOS", "1.0.0", ids[1],
		"text/xml;subtype=&quot;om/1.0.0&quot;")
getobsbyid
encodeRequestXML(getobsbyid) # is valid!

getobsbyid <- GetObservationById(service = "SOS", version = "1.0.0",
		observationId = ids[2],
		responseFormat = "text/xml;subtype=&quot;om/1.0.0&quot;",
		resultModel = "om:Measurement", responseMode = "inline",
		srsName = "urn:ogc:def:crs:EPSG:4326")
getobsbyid
encodeRequestXML(getobsbyid) # is valid!

weathersos = SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos",
		method = "POST", verboseOutput = TRUE)
ids = paste("o_33765", c(80:99), sep = "")
obs <- getObservationById(sos = weathersos, observationId = ids[[1]],
		resultModel = SosSupportedResultModels()[[2]])
str(obs@result)
# om:Observation works!

meas <- getObservationById(sos = weathersos, observationId = ids[[1]],
		resultModel = SosSupportedResultModels()[[1]])
str(meas@result)
# om:Measurement works!

obs <- getObservationById(sos = weathersos, observationId = ids[[1]], verbose = TRUE)
obs <- getObservationById(sos = weathersos, observationId = ids[[1]], resultModel = SosSupportedResultModels()[[2]], verbose = TRUE)
# works!

################################################################################
# temporal operations
weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos", method = "POST")

#format(as.POSIXct("2010-07-01 12:00"), weathersos@timeFormat)
t1 <- sosCreateTimeInstant(sos = weathersos, time = as.POSIXct("2010-07-01 12:00"), frame = "frameLala", calendarEraName = "calName", indeterminatePosition = "YES")
p1 <- sosCreateTimePeriod(sos = weathersos, begin = as.POSIXct("2010-03-01 12:15"),
		end = as.POSIXct("2010-03-02 12:15"), timeInterval = GmlTimeInterval(interval = "1d", unit = "day",
				radix = as.integer(17), factor = as.integer(42)))
# works!

eventTime1 <- sosCreateEventTimeList(op = SosSupportedTemporalOperators()[["TM_During"]], t = p1)
eventTime2 <- sosCreateEventTimeList(op = SosSupportedTemporalOperators()[["TM_Equals"]], t = t1)
eventTime3 <- sosCreateEventTimeList(op = SosSupportedTemporalOperators()[["TM_After"]], t = t1)
eventTime4 <- sosCreateEventTimeList(op = SosSupportedTemporalOperators()[["TM_Before"]], t = t1)
encodeXML(eventTime3[[1]], sos = weathersos, v = TRUE)
encodeXML(eventTime4[[1]], sos = weathersos)
encodeXML(eventTime2[[1]], sos = weathersos)
encodeXML(eventTime1[[1]], sos = weathersos, v = TRUE)
encodeKVP(eventTime1[[1]], sos = weathersos)
encodeKVP(eventTime2[[1]], sos = weathersos)
# works!

################################################################################
source("/home/daniel/Dokumente/2010_SOS4R/workspace/sos4R/sandbox/loadSources.R")
################################################################################

################################################################################
# spatial operations

# valid foi ids for WeatherSOS:
foiIDs <- list(
		"urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93",
		"urn:ogc:object:feature:OSIRIS-HWS:efeb807b-bd24-4128-a920-f6729bcdd111",
		"id:02")
foiObj <- sosCreateFeatureOfInterest(objectIDs = foiIDs[1:2])
encodeXML(foiObj)
# ok

bbox1 <- sosCreateBBoxMatrix(lowLat = 50.0, lowLon = 7.0, uppLat = 53.0, uppLon = 10.0)
foiBBox1 <- sosCreateFeatureOfInterest(bbox = bbox1, srsName = "urn:ogc:def:crs:EPSG:6.8:4326")
encodeXML(foiBBox1, verbose = TRUE)
# ok

foiBBox2 <- sosCreateFeatureOfInterest(
		spatialOps = sosCreateBBOX(lowLat = 50.0, lowLon = 7.0,
				uppLat = 53.0, uppLon = 10.0, 
				srsName = "urn:ogc:def:crs:EPSG:6.8:4326",
				srsDimension = as.integer(2), axisLabels = "lat,lon",
				uomLabels = "deg,deg",
				propertyName = "lalaUndPooh"))
encodeXML(foiBBox2, weathersos)
# ok!

################################################################################
source("/home/daniel/Dokumente/2010_SOS4R/workspace/sos4R/sandbox/loadSources.R")
################################################################################

################################################################################
# GetObservations
weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")

# LATEST
obs <- getObservation(sos = weathersos,
		observedProperty = sosObservedProperties(weathersos)[4],
		offering = sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]],
		latest = TRUE,
		inspect = TRUE, verbose = TRUE)
obs@result
# works!

# TIME INTERVAL
go.eventTime1 = sosCreateEventTimeList(
		operator = SosSupportedTemporalOperators()[["TM_After"]],
		sosCreateTimeInstant(sos = weathersos,
				time = as.POSIXct("2010-09-20 18:00")))
# TM_After does not seem to work, probably error in sos?!

go.eventTime1a = sosCreateEventTimeList(
		operator = SosSupportedTemporalOperators()[["TM_Equals"]],
		sosCreateTimeInstant(sos = weathersos,
				time = as.POSIXct("2010-09-20 18:00")))
go.eventTime1b = sosCreateEventTimeList(
		operator = SosSupportedTemporalOperators()[["TM_Equals"]],
		sosCreateTimeInstant(sos = weathersos,
				time = as.POSIXct("2010-09-20 18:15")))
obs2 <- getObservation(sos = weathersos,
		observedProperty = list(go.observedProperty),
		procedure = list(sosProcedures(weathersos)[[1]]),
		eventTime = list(go.eventTime1a, go.eventTime1b),
		offering = go.offering@id, inspect = TRUE)
obs2[[1]]@result
# weird, as the sos returns two equal observations here...

go.eventTime3 = sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
				begin = as.POSIXct("2010-09-16 18:00"),
				end = as.POSIXct("2010-09-20 18:00")))
obs3 <- getObservation(sos = weathersos,
		observedProperty = list(go.observedProperty),
		procedure = list(sosProcedures(weathersos)[[1]]),
		eventTime = go.eventTime3,
		offering = go.offering@id)
#obs3@result
# heureka!
summary(obs3@result) # finally!
plot(x = obs3@result[["Time"]],
		y = obs3@result[["urn:ogc:def:property:OGC::Temperature"]],
		type = "l",
		main = "Temperature in Muenster",
		xlab = "Time",
		ylab = "Temperature (degree C)")


# two procedures
obs4 <- getObservation(sos = weathersos,
		observedProperty = sosObservedProperties(weathersos)[5],
		procedure = sosProcedures(weathersos),
		eventTime = sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
						begin = as.POSIXct("2009-08-10 12:00"),
						end = as.POSIXct("2009-08-20 12:00"))),
		#featureOfInterest = foiBBox,
		offering = sosOfferings(weathersos)[[5]])
str(obs4[[1]]@result)
str(obs4[[2]]@result)
# Can I automatically join these? No, not really, as time stamps differ!
data.frame(obs4[[1]]@result["Time"][1:10,], obs4[[2]]@result["Time"][1:10,])

# Attention: plots ignore the fact that the times do NOT perfectly match!
x <- 800
plot(x = obs4[[1]]@result[[1]][1:x], y = obs4[[1]]@result[[3]][1:x], type = "l",
		col = "steelblue", main = "Temperature in Muenster and Kaernten, 2009",
		xlab = "Time (00:00 o'clock)",
		ylab = "Temperature (degree C)",
		xaxt="n") # do not ploplott x-axis
r <- as.POSIXct(round(range(obs4[[1]]@result[[1]]), "days"))
axis.POSIXct(side = 1, x = obs4[[1]]@result[[1]][1:x], format = "%d. %h",
		at = seq(r[1], r[2], by="day"))
lines(x = obs4[[2]]@result[[1]][1:x], y = obs4[[2]]@result[[3]][1:x],
		col = "orange")
legend("topleft", legend = c("Muenster", "Kaernten"),
		col = c("steelblue", "orange"), lty = 1, bty="n")

savePlot(type = "png", filename = "/tmp/temp-MS-K.png")

hist(obs4[[2]]@result[[3]])
pairs(obs4[[2]]@result) # not really useful :-)

# FOI with ids, varying @id and if procedures are given or not
foiIDs <- sosCreateFeatureOfInterest(
		objectIDs = list(
				"urn:ogc:object:feature:OSIRIS-HWS:3d3b239f-7696-4864-9d07-15447eae2b93",
				"urn:ogc:object:feature:OSIRIS-HWS:efeb807b-bd24-4128-a920-f6729bcdd111",
				"id:02")[1:1])
obs5 <- getObservation(sos = weathersos,
		observedProperty = sosObservedProperties(weathersos)[2],
		#procedure = sosProcedures(weathersos),
		eventTime = sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
						begin = as.POSIXct("2009-08-10 12:00"),
						end = as.POSIXct("2009-08-12 12:00"))),
		featureOfInterest = foiIDs,
		offering = sosOfferings(weathersos)[[2]], # @id not needed, works!
		inspect = TRUE)
str(obs5)
# only returns the (one) selected feature
# works!

# BBOX, different variants, also testing more several phenomena, but in
# WeatherSOS, there is only one observedProperty per offering
foiBBoxMS <- sosCreateFeatureOfInterest(
		bbox = sosCreateBBoxMatrix(lowLat = 50.0, lowLon = 7.0,
				uppLat = 53.0, uppLon = 10.0),
		srsName = "urn:ogc:def:crs:EPSG:4326")
foiBBoxMSandK <- sosCreateFeatureOfInterest(
		bbox = sosCreateBBoxMatrix(lowLat = 45.0, lowLon = 7.0,
				uppLat = 53.0, uppLon = 15.0),
		srsName = "urn:ogc:def:crs:EPSG:4326")
foiBBoxNone <- sosCreateFeatureOfInterest(
		bbox = sosCreateBBoxMatrix(lowLat = 0.0, lowLon = 0.0,
				uppLat = 1.0, uppLon = 1.0),
		srsName = "urn:ogc:def:crs:EPSG:4326")
obs6 <- getObservation(sos = weathersos,
		observedProperty = sosObservedProperties(weathersos)[2:5],
		procedure = sosProcedures(weathersos),
		eventTime = sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
						begin = as.POSIXct("2009-08-10 12:00"),
						end = as.POSIXct("2009-08-17 12:00"))),
		featureOfInterest = foiBBoxMSandK,
		offering = sosOfferings(weathersos)[[3]]@id, inspect = TRUE,
		resultModel = SosSupportedResultModels()[[1]])
str(obs6)
# works!

obs6@result[[1]][c(1,672)]; obs6@samplingTime#  both the same values - good!

# testing with om:Measurement
meas1 <- getObservation(sos = weathersos,
		observedProperty = sosObservedProperties(weathersos),
		procedure = sosProcedures(weathersos),
		eventTime = sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
						begin = as.POSIXct("2009-08-10 12:00"),
						end = as.POSIXct("2009-08-11 12:00"))),
		offering = sosOfferings(weathersos)[[3]]@id, inspect = TRUE,
		resultModel = SosSupportedResultModels()[[1]])
# get values from measurements
values1 <- sapply(sapply(meas1, slot, name = "result"), slot, name = "value")
values1 <- as.double(values1)
times1 <- lapply(sapply(sapply(meas1, slot, name = "samplingTime"), 
				slot, name = "timePosition"), slot, name = "time")

df.times1 <- as.data.frame(times1); names(df.times1) <- "time"
df.values1 <- as.data.frame(values1); names(df.values1) <- "values"

data.frame(times1, values1)
plot(x = times1, y = values1)

# two observedPropertys is not possible with weathersos, as the offerings contain only one phenomenon each

################################################################################
# ogc:comparisonOps dummy classes and workaround for passing though XML
manualResult <- '<sos:result><ogc:PropertyIsGreaterThan><ogc:PropertyName>urn:ogc:def:phenomenon:OGC:1.0.30:waterlevel</ogc:PropertyName><ogc:Literal>5</ogc:Literal></ogc:PropertyIsGreaterThan></sos:result>'
manualResult1 <- xmlParseString(manualResult)
class(manualResult1)
# [1] "XMLInternalDocument" "XMLAbstractDocument" "oldClass"  

manualResult2 <- xmlNode(name = sosResultName)
manualResult2$children[[1]] <- xmlNode(name = ogcComparisonOpEqualToName)
class(manualResult2)
# [1] "XMLNode" "RXMLAbstractNode" "XMLAbstractNode"  "oldClass"

manualResult3 <- xmlParse(manualResult)
class(manualResult3)
# [1] "XMLInternalDocument" "XMLAbstractDocument" "oldClass"

encodeXML(manualResult1, verbose = TRUE)
encodeXML(manualResult2, verbose = TRUE)
encodeXML(manualResult3, verbose = TRUE)
# works, warnings as well.

weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
# this now even works with defaults for procedures and observerProperty, i.e.
# just taking all available!
lastTenHours <- sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
		begin = (Sys.time() - 36000), end = Sys.time()))
obs1 <- getObservation(sos = weathersos, eventTime = lastTenHours,
		offering = sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]])

myResult1 <- xmlParseString('<sos:result><ogc:PropertyIsGreaterThan>
<ogc:PropertyName>urn:ogc:def:property:OGC::Temperature</ogc:PropertyName>
<ogc:Literal>20</ogc:Literal></ogc:PropertyIsGreaterThan></sos:result>')
# results in problems with external pointers...

pn <- xmlNode(name = ogcPropertyNameName, namespace = ogcNamespacePrefix)
xmlValue(pn) <- "urn:ogc:def:property:OGC::Temperature"
l <- xmlNode(name = "Literal", namespace = ogcNamespacePrefix)
xmlValue(l) <- "10"
comp <- xmlNode(name = ogcComparisonOpGreaterThanName,
		namespace = ogcNamespacePrefix, .children = list(pn, l))
myResult2 <- xmlNode(name = sosResultName, namespace = sosNamespacePrefix,
		.children = list(comp))

obs2 <- getObservation(sos = weathersos, eventTime = lastTenHours,
		offering = sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]],
		result = myResult2, verbose = TRUE)

length(sosResult(obs1)[[3]]); length(sosResult(obs2)[[3]])
# WORKS!


################################################################################
# using attributes on data frame
tempDataFrame <- sosResult(pegelObs)
attributes(tempDataFrame)
#$names
#[1] "Time"             "feature"          "Wasserstand [cm]"
#$class
#[1] "data.frame"
#$row.names
#[1]    1    2    3    4    5    6    7    8    9   10   11   12   13   14

attributes(tempDataFrame[,1])
#$class
#[1] "POSIXt"  "POSIXct"
#$tzone
#[1] ""

attributes(tempDataFrame[,2])
#$levels
#[1] "Steinbach_2450010"
#$class
#[1] "factor"

attributes(tempDataFrame[,3])
#NULL

tempDataFrame[1:2,]
#				 Time           feature Wasserstand [cm]
#1 2010-08-22 00:15:00 Steinbach_2450010              171
#2 2010-08-22 00:30:00 Steinbach_2450010              170


################################################################################
# Was I stupid?
weathersos = SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
sensor <- describeSensor(weathersos, sosProcedures(weathersos)[[1]])
sensorMLtemplate <- makeClassTemplate(xnode = xmlRoot(sensor@xml))
str(sensorMLtemplate)
cat(sensorMLtemplate$def)
cat(sensorMLtemplate$slots[["member"]])
cat(sensorMLtemplate$coerce)

################################################################################
# Added OmObservationCollection
pegelsos <- SOS(url = "http://v-sos.uni-muenster.de:8080/PegelOnlineSOSv2/sos")
pegelObs <- getObservation(sos = pegelsos,
		observedProperty = sosObservedProperties(pegelsos)[3],
		offering = sosOfferings(pegelsos)[[1]],
		procedure = sosProcedures(pegelsos)[c(2501,2503,2505)],
		eventTime = sosCreateEventTimeList(time = sosCreateTimePeriod(
						sos = pegelsos,
						begin = Sys.time() - (3600 * 24), # * 360),
						end = Sys.time())))

sosResult(pegelObs)[[1]][1,]
# is the same as
sosResult(pegelObs@members[[1]])[1,]

length(pegelObs)

# added indexing functions
pegelObs[1]
str(pegelObs[1])
str(pegelObs[1:2])
str(pegelObs[c(3,4)])
str(pegelObs[c(3:6)])


################################################################################
# coercion
str(as.list(pegelObs)[3:5])
str(as.data.frame(pegelObs[[1]]))
# works!

################################################################################
# sosFeaturesOfInterest and sosOfferings
str(sosFeaturesOfInterest(sosOfferings(weathersos)[[1]]))
sosOfferings(weathersos, c("WIND_DIRECTION", "ATMOSPHERIC_PRESSURE"))
sosOfferingIds(weathersos)

tempIds <- sosOfferingIds(weathersos)
tempOfferings <- sosOfferings(weathersos)
sosObservedProperties(sosOfferings(weathersos)[[1]])
sosObservedProperties(weathersos)

################################################################################
# encode xml character string (again)
isXMLString('<lala name="horst"><pooh /></lala>')
xmlParseString('<lala name="horst"><pooh /></lala>')
# ok... put into encodeXML

xml1 <- encodeXML('<lala name="horst"><pooh /></lala>', weathersos)
str(xml1)
# good

pn <- xmlNode(name = ogcPropertyNameName, namespace = ogcNamespacePrefix)
xmlValue(pn) <- "urn:ogc:def:property:OGC::Temperature"
l <- xmlNode(name = "Literal", namespace = ogcNamespacePrefix)
xmlValue(l) <- "3"
comp <- xmlNode(name = ogcComparisonOpLessThanOrEqualToName,
		namespace = ogcNamespacePrefix, .children = list(pn, l))
myResult2 <- xmlNode(name = sosResultName, namespace = sosNamespacePrefix,
		.children = list(comp))
str(myResult2)
class(myResult2)
lastTenHours <- sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
				begin = (Sys.time() - 36000), end = Sys.time()))
obs2 <- getObservation(sos = weathersos, eventTime = lastTenHours,
		offering = sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]],
		result = myResult2, inspect = TRUE)
# works!


sosResultString <- '<sos:result><ogc:PropertyIsGreaterThan>
		<ogc:PropertyName>urn:ogc:def:property:OGC::Temperature</ogc:PropertyName>
		<ogc:Literal>20</ogc:Literal></ogc:PropertyIsGreaterThan></sos:result>'
myResult1 <- xmlParseString(sosResultString, clean = FALSE)
# namespace prefixes not found

shouldWorkString <- '<sos:result xmlns:sos="http://www.opengis.net/sos/1.0"
		xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance"
		xmlns:ows="http://www.opengis.net/ows/1.1"
		xmlns:om="http://www.opengis.net/om/1.0"
		xmlns:ogc="http://www.opengis.net/ogc"
		xmlns:gml="http://www.opengis.net/gml">
		<ogc:PropertyIsGreaterThan>
		<ogc:PropertyName>urn:ogc:def:property:OGC::Temperature</ogc:PropertyName>
		<ogc:Literal>20</ogc:Literal></ogc:PropertyIsGreaterThan></sos:result>'
xmlParseString(shouldWorkString)
str(encodeXML(shouldWorkString, weathersos))

# try automatic namespace adding with encodeXML
namespacedResult <- xmlNode(name = sosResultName,
		namespace = sosNamespacePrefix,
		namespaceDefinitions = c(.sosNamespaceDefinitionsForAll,
				.sosNamespaceDefinitionsGetObs))
namespacedString <- sub(pattern = result, replacement = paste("result"))
# all good.

# some trick with handlers  could work ...
xmlParse(sosResultString, asText = TRUE, handlers=list("result"=function(x, ...){xmlParse(x)}))
# but doesn't right now

.hack <- 'result xmlns:sos="http://www.opengis.net/sos/1.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ows="http://www.opengis.net/ows/1.1" xmlns:om="http://www.opengis.net/om/1.0" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml">'
hackedString <- sub(pattern = "result", replacement = .hack, x = sosResultString)
xmlParseString(hackedString, clean = FALSE)
# all better

# move hack to encodeXML
encodeXML(obj = sosResultString, sos = weathersos)
.result <- encodeXML(obj = sosResultString, sos = weathersos, addNamespaces = TRUE)
# cannot attach .result this to anything!
n <- xmlNode("lala")
n[[1]] <- list(.result)
n <- addChildren(node = n, kids = list(.result))
#Fehler in as.vector(x, "list") : 
#		cannot coerce type 'externalptr' to vector of type 'list'

lastTenHours <- sosCreateEventTimeList(sosCreateTimePeriod(sos = weathersos,
				begin = (Sys.time() - 36000), end = Sys.time()))
obs2 <- getObservation(sos = weathersos, eventTime = lastTenHours,
		offering = sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]],
		result = sosResultString, inspect = TRUE)
# STILL extrnlpntr problem
debug(.sosEncodeRequestXMLGetObservation_1.0.0)

obs2 <- getObservation(sos = weathersos, eventTime = lastTenHours,
		offering = sosOfferings(weathersos)[["ATMOSPHERIC_TEMPERATURE"]],
		result = xmlParseString(hackedString), inspect = TRUE)

################################################################################
# accessing feature and it's coordinates
obs <- pegelObs[[1]]
obs@featureOfInterest@featureMembers[[1]]@feature@position@point@pos@pos

str(sosFeaturesOfInterest(obs))

sosCoordinates(obs@featureOfInterest@featureMembers[[1]]@feature@position@point@pos)
sosCoordinates(obs@featureOfInterest@featureMembers[[1]]@feature@position@point)
sosCoordinates(obs@featureOfInterest@featureMembers[[1]]@feature@position)
sosCoordinates(obs@featureOfInterest@featureMembers[[1]]@feature)
sosCoordinates(obs@featureOfInterest@featureMembers[[1]])
sosCoordinates(obs@featureOfInterest@featureMembers)
sosCoordinates(obs@featureOfInterest)
sosCoordinates(obs)

sosCoordinates(pegelObs)
sosCoordinates(pegelObs[[1]])
sosCoordinates(pegelObs[[3]])

sosCoordinates(observation.pm10.week[1:2])
sosCoordinates(observation.pm10.week)
# good so far

# add the feature, id is in sampling point
foiid <- str(observation.pm10.week[[1]]@featureOfInterest@featureMembers[[1]]@feature@id)
foiid
sosId(observation.pm10.week[[1]]@featureOfInterest@featureMembers[[1]]@feature)

#
sosSrsName(obs@featureOfInterest@featureMembers[[1]]@feature@position@point@pos)
sosSrsName(obs@featureOfInterest@featureMembers[[1]]@feature@position@point)

################################################################################
# indexing OmObservationCollection with observed properties and procedures
.getObservationsWithProcedure(observation.pm10.week,
		sosProcedures(observation.pm10.week)[[5]])
.getObservationsWithProcedure(observation.pm10.week,
		sosProcedures(observation.pm10.week)[2:4])

observation.pm10.week[which(sosProcedures(observation.pm10.week) %in% 
				c("urn:x-eea:object:sensor:airquality:PL0014A",
						"urn:x-eea:object:sensor:airquality:AT30599"))]
observation.pm10.week[which(sosProcedures(observation.pm10.week) %in% 
						"urn:x-eea:object:sensor:airquality:AT30599")]

length(.getObservationsWithObservedProperty(observation.pm10.week,
		sosObservedProperties(observation.pm10.week)))

observation.pm10.week["urn:x-eea:object:sensor:airquality:CZ0HHKB"]
observation.pm10.week[sosProcedures(observation.pm10.week)[4:6]]
observation.pm10.week[sosObservedProperties(observation.pm10.week)[1]]
# works!

.getObservationsWithFoiId(observation.pm10.week, "foi_CZ0UDCM")

observation.pm10.week[sosProcedures(observation.pm10.week)[1]]
observation.pm10.week[sosProcedures(observation.pm10.week)[[1]]]
observation.pm10.week["foi_CZ0UDCM"]

################################################################################
# binding result data frames from a observation collection
sosResult(observation.pm10.week[1:3])

# try that with different observed properties
pegelsos <- SOS(url = "http://v-sos.uni-muenster.de:8080/PegelOnlineSOSv2/sos")

sosOfferings(pegelsos)
offering.wasserstand <- sosOfferings(pegelsos)[["WASSERSTAND_ROHDATEN"]]

lastHour <- sosCreateEventTimeList(sosCreateTimePeriod(sos = pegelsos,
				begin = (Sys.time() - 3600), end = Sys.time()))
pegelObs.10 <- getObservation(sos = pegelsos,
		observedProperty = sosObservedProperties(offering.wasserstand)[5],
		offering = wasserstand_roh,
		procedure = c("Wasserstand-Bake_Z_9510066",
				"Wasserstand-Bake_A_9510063"),
		eventTime = lastHour,
		verbose = TRUE)
# this does not return the requested observed property, but just wasserstand...


myOff <- sosOfferings(MBARI)[[1]]
myProc <- sosProcedures(MBARI)[[1]]
mbariObs <- getObservation(sos = MBARI, offering = myOff, procedure = myProc)
sosResult(mbariObs[[1]])[1:2,]
one <- sosResult(mbariObs[[1]])
two <- one[,1:5]
do.call(rbind, list(one, two))[1:2,]
# does NOT work with different number of columns!

################################################################################
# plotting offerings
library(sp)

weathersos <- SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
weathersos.offerings <- sosOfferings(weathersos)

off.ids <- sosId(weathersos.offerings)
off.bounds <- sosBoundedBy(weathersos.offerings)

ll <- as.numeric(strsplit(off.bounds[[1]]$lowerCorner, " ")[[1]])
uu <- as.numeric(strsplit(off.bounds[[1]]$upperCorner, " ")[[1]])
ll; uu

off.1.poly <- cbind(c(ll[[2]], ll[[2]], uu[[2]], uu[[2]], ll[[2]]),
				c(ll[[1]], uu[[1]], uu[[1]], ll[[1]], ll[[1]]))
off.1.spatPoly <- SpatialPolygons(list(Polygons(
						list(Polygon(off.1.poly)), off.ids[[1]])))

off.1.epsg <- strsplit(off.bounds[[1]]$srsName, ":")[[1]]
off.1.epsg <- off.1.epsg[[length(off.1.epsg)]]

proj4string(off.1.spatPoly) <- paste("+init=epsg:", off.1.epsg, sep = "")

off.1.spdf <- SpatialPolygonsDataFrame(off.1.spatPoly,
		as.data.frame(c(42), row.names = c(off.ids[[1]])))

off.1.spatPoly
off.1.spdf

plot(off.1.spdf)
# plots bbox!


################################################################################
# saving of the orginal document
weathersos <- SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
weathersos.offerings <- sosOfferings(weathersos)

# by id
obs <- getObservationById(sos = weathersos, observationId = "o_3508493",
		saveOriginal = "lalalaaaaa", verbose = TRUE)

.begin <- sosTime(weathersos.offerings[[1]], convert = TRUE)[["begin"]]
.timelist <- sosCreateEventTimeList(sosCreateTimePeriod(weathersos,
				begin = .begin,
				end = (.begin + 3600*24)))
obs <- getObservation(sos = weathersos, offering = weathersos.offerings[[1]],
		eventTime = .timelist, saveOriginal = TRUE, verbose = TRUE)
obs <- getObservation(sos = weathersos, offering = weathersos.offerings[[1]],
		eventTime = .timelist, saveOriginal = "text", verbose = TRUE)

# read it back
.filename <- paste("o_3508493", , ".xml", sep = "")
obs.2 <- sosParse(weathersos, .filename)

describeSensor(sos = weathersos,
		procedure = sosProcedures(weathersos)[[1]][[1]],
		saveOriginal = TRUE)
describeSensor(sos = weathersos,
		procedure = sosProcedures(weathersos)[[1]][[1]],
		saveOriginal = "sensor11.xml")

################################################################################
# url from DCPs
weathersos <- SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
weathersos.offerings <- sosOfferings(weathersos)

# manually
weathersos@capabilities@operations@operations$GetCapabilities@DCPs["Get"]

sosGetDCP(weathersos, operation = sosGetObservationName)
sosGetDCP(weathersos, operation = sosGetObservationName, type = "Get")
sosGetDCP(weathersos, operation = sosGetObservationName, type = "Post")

describeSensor(weathersos, sosProcedures(weathersos)[[1]][[1]])
getObservation(weathersos, sosOfferings(weathersos)[[1]], latest = TRUE)
# still works

################################################################################
# plotting and coercion
csiro <- SOS("http://wron.net.au/CSIRO_SOS/sos")
rainfall.off.csiro <- sosOfferings(csiro)["Rain Gauges"][[1]]
sosBoundedBy(rainfall.off.csiro, bbox = TRUE)
sosCapabilitiesDocumentOriginal(csiro)
# wrong coordinate order in capabilities!!


csiro <- SOS("http://wron.net.au/CSIRO_SOS/sos", switchCoordinates = TRUE)
rainfall.off.csiro <- sosOfferings(csiro)["Rain Gauges"][[1]]
sosBoundedBy(rainfall.off.csiro, bbox = TRUE)
plot(rainfall.off.csiro)

poly <- as(rainfall.off.csiro, "SpatialPolygons")
poly
poly2 <- as(rainfall.off.bom, "SpatialPolygons")
bbox(poly)
map.where("world", coordinates(poly))

# get the map data dir
Sys.getenv("R_MAP_DATA_DIR")
# check out world.N

weathersos <- SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
weathersos.offerings <- sosOfferings(weathersos)

plot(weathersos.offerings[[1]])
plot(weathersos, regions = c("Germany", "Austria"))
plot(csiro)

sosBoundedBy(weathersos.offerings[[1]])
sosBoundedBy(weathersos.offerings[[1]], bbox = TRUE)
as(weathersos.offerings[[1]], "SpatialPolygons")

#
# plotting sensor positions
#
proc1 <- sosProcedures(weathersos)[[1]][[1]]
proc1.descr <- describeSensor(weathersos, proc1, verbose = TRUE)

# weathersos
sosId(proc1.descr)
sosName(proc1.descr) # short name
sosAbstract(proc1.descr) # description name
coords <- sosCoordinates(proc1.descr, sos = weathersos, verbose = TRUE)
coords
attributes(coords[,1])
attributes(coords[,3])

proc.all <- sosProcedures(weathersos)[[1]]
proc.all.descr <- lapply(proc.all, describeSensor, sos = weathersos)
coords.all <- sosCoordinates(proc.all.descr, sos = weathersos)
coords.all
attributes(coords.all[,3])

# try csiro, but not successful
proc2 <- sosProcedures(csiro)[[3]][[2]]
proc2.descr <- describeSensor(sos = csiro, procedure = proc2)
proc2.descr@xml
sosId(proc2.descr)
sosName(proc2.descr)
sosAbstract(proc2.descr)

###
# TODO continue here with plotting of sensor
weathersos <- SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
proc1 <- sosProcedures(weathersos)[[1]][[1]]
proc1.descr <- describeSensor(weathersos, proc1, verbose = TRUE)

# convert to spatial points data frame
as(proc1.descr, "Spatial")

.coords <- attributes(proc1.descr@coords)
.crds <- .coords[,c("x", "y")]
.data <- .coords[,!colnames(.coords)%in%c("x", "y")]
.crs <- attributes(.coords)[["referenceFrame"]]

SpatialPointsDataFrame(coords = .crds, data = .data, crs = .crs)

summary(coords)

# this works, move it to plot function...
plot(coords, add = TRUE, pch = 19)
text(x = coordinates(coords)[,"x"], y = coordinates(coords)[,"y"],
		labels = row.names(coords@data), adj = c(0, 1), cex = 0.75)


################################################################################
# plot sensor positions -- SEE WEATHERSOS DEMO


################################################################################
# fixing errors in vignette...
weathersos <- SOS(url = "http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")

sosFeaturesOfInterest(weathersos)[1:2]

################################################################################
# requesting offerings by name
aqe.converters <- SosDataFieldConvertingFunctions(
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[PM10]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[NO2]" = sosConvertDouble,
		"http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[O3]" = sosConvertDouble,
		"http://www.opengis.net/def/property/OGC/0/SamplingTime" = sosConvertTime,
		"http://www.opengis.net/def/property/OGC/0/FeatureOfInterest" = sosConvertString)
aqe <- SOS(url = "http://giv-uw.uni-muenster.de:8080/AQE/sos",
		dataFieldConverters = aqe.converters)
no2off <- sosOfferings(aqe)[["NO2"]]

prop <- "http://giv-genesis.uni-muenster.de:8080/SOR/REST/phenomenon/OGC/Concentration[NO2]"
#foi = "foi_DEST080"
foi <- sosCreateFeatureOfInterest(sosFeaturesOfInterest(no2off)[1:20])
time <- sosCreateEventTimeList(sosCreateTimePeriod(sos = aqe,
				begin = as.POSIXct("2008-12-21T23:01:00Z"),
				end = as.POSIXct("2008-12-31T23:01:00Z")))

obs <- getObservation(sos = aqe, offering = "NO2", verbose = TRUE,
#		observedProperty = sosObservedProperties(sosOfferings(aqe)[["NO2"]]),
		featureOfInterest = foi,
		eventTime = time)


################################################################################
# r-help with time zones

sessionInfo()

t1 <- strptime("1995-05-25T15:30:00+10:00", format = "%Y-%m-%dT%H:%M:%OS")
t2 <- strptime("1995-05-25T15:30:00+10:00", format = "%Y-%m-%dT%H:%M:%OS%z")

strftime(t1, format = "%Y-%m-%dT%H:%M:%OS")
strftime(t1, format = "%Y-%m-%dT%H:%M:%OS%z")
# Ends in "Mitteleuropaeische Sommerzeit", not in +10:00, so time zone is ignored!
# Also no difference beetween %z and %z !
strftime(t1, format = "%Y-%m-%dT%H:%M:%OS%Z")
# All this does NOT remove the "Mitteleuropaeische Zeit" from the strftime output!!

# Can locale solve the problem?
Sys.getlocale(category = "LC_TIME")
Sys.setlocale("LC_TIME", "English")

strftime(t1, format = "%Y-%m-%dT%H:%M:%OS%z")
# [1] "1995-05-25T15:30:00Mitteleuropaeische Sommerzeit" -- No change.

# does t1 actually have time zone?
attributes(t1)

format(t1, format = "%Y-%m-%dT%H:%M:%OS%z") # usetz = TRUE) # no change on usetz

# Is the : in offset the problem?
t3 <- strptime("1995-05-25T15:30:00+1000", format = "%Y-%m-%dT%H:%M:%S%z")
attributes(t3)
format(t3, format = "%Y-%m-%dT%H:%M:%OS%z")
# [1] "1995-05-25T07:30:00Mitteleuropaeische Sommerzeit"

strftime(t1, format = "%Y-%m-%dT%H:%M:%OS%z", tz = "+0200") # no effect on setting tz

Sys.setenv(TZ="GMT") # no working effect on format and strftime


###############
# http://r.789695.n4.nabble.com/Timezone-issue-with-strftime-strptime-and-z-and-Z-tt3346204.html#none

# SECOND EMAIL DAVID
as.POSIXlt(gsub("T", " ", 	# change T to space
							# but preserve the sign for the %z format string
				gsub("(T..:..:.....):", "\\1", "1995-05-25T15:30:00-10:00")),
		format="%Y-%m-%d %H:%M:%S%z")
# [1] "1995-05-25 21:30:00"

# To get output in GMT add tz argument to as.POSIXlt:
as.POSIXlt(gsub("T", " ", 	# change T to space
							# but preserve the sign for the %z format string
				gsub("(T..:..:.....):", "\\1", "1995-05-25T15:30:00-10:00")),
		format="%Y-%m-%d %H:%M:%S%z", tz="GMT")

# Does still not ouput the numerical time zones!

# To get output in GMT add tz argument to as.POSIXlt:
time <- as.POSIXlt(gsub("T", " ", 	# change T to space
				# but preserve the sign for the %z format string
				gsub("(T..:..:.....):", "\\1", "1995-05-25T15:30:00-10:00")),
		format="%Y-%m-%d %H:%M:%S%z") #, tz="GMT")
format(x = time, format = "%Y-%m-%d %H:%M:%S%z")
# [1] "1995-05-26 01:30:00Mitteleuropaeische Zeit"
strftime(x = time, format = "%Y-%m-%d %H:%M:%S%z")
# [1] "1995-05-26 03:30:00Mitteleuropaeische Sommerzeit"

x <- as.POSIXlt(gsub("T", " ", #change T to space
				# but preserve the sign for the %z format string
		gsub("(T..:..:.....):", "\\1", "1995-05-25T15:30:00-10:00")),
		format="%Y-%m-%d %H:%M:%S%z", tz="GMT")
x
# [1] "1995-05-26 01:30:00 GMT"
format(x, "%Y-%m-%d %H:%M:%S%z")
# [1] "1995-05-26 03:30:00Mitteleuropaeische Zeit"
format(x, "%Y-%m-%dT%H:%M:%S%z")
# [1] "1995-05-26T01:30:00Mitteleuropaeische Zeit"


################################################################################
#
# sos4R test script EDC by EP, to be run after edc-entwicklerforum.R
library(gstat)

# move to an appropriate CRS:
RD = CRS(paste("+init=epsg:28992",
				"+towgs84=565.237,50.0087,465.658,-0.406857,0.350733,-1.87035,4.0812"))

library(rgdal)
no2_spdf = spTransform(no2_spdf, RD)
map.lines = spTransform(map.lines, RD)

no2.T1 = no2_spdf[no2_spdf$SamplingTime == min(no2_spdf$SamplingTime),]

grd = SpatialPixels(SpatialPoints(makegrid(bbox(map.lines), n = 1000)),
		proj4string = proj4string(map.lines))

plot(grd)
plot(map.lines, add=T)

names(no2.T1)[3] = "NO2"
NO2.idw =idw(NO2~1, no2.T1, grd)
lt = list(list("sp.lines", map.lines),
		list("sp.points", no2.T1, col = grey(.5)))
spplot(NO2.idw[1], col.regions = bpy.colors(), sp.layout=lt)

# plot stations, something wrong with projection?
plot(no2_spdf, add = TRUE)


################################################################################
# Event time list creation based on character strings:

sosCreateTime(sos = aqe, time = "2007-08-01 08:00::2007-08-05 15:00")
sosCreateTime(sos = aqe, time = "2007-08-01 15:00/2007-08-05 20:00")

# test with request
getObservation(sos = aqe, offering = sosOfferings(aqe)[[1]],
		eventTime = sosCreateTime(sos = aqe,
				time = "2007-08-01 15:00/2007-08-05 20:00"))
# all good!

sosCreateTime(sos = aqe, time = "::2007-08-05")
sosCreateTime(sos = aqe, time = "2007-08-05::")

getObservation(sos = aqe, offering = sosOfferings(aqe)[[1]],
		eventTime = sosCreateTime(sos = aqe, time =  "::2007-08-02"))

sosCreateTime(sos = aqe, time = "/2007-08-05")
sosCreateTime(sos = aqe, time = "2007-08-05/")


################################################################################
# get all code chunks of a vignette
edit(vignette("sos4R"))


################################################################################
# there are some characters that should not be used for colum names as they can
# cause problems when using formula (?formula)
.escapeColumnName("Concentration[NO2]")
.escapeColumnName("Concentration[[NO2][und@home.net]]")
.escapeColumnName("Concentration~A+B-C")
.escapeColumnName("Concentration**A$A$A**")
# works.

################################################################################
# new method sosGetUOM

sosUOM(GmlMeasure(42.0, "m"))

sosUOM(obs.temp.latest[[1]])
sosUOM(obs.temp.latest[1:2])
sosUOM(obs.temp.latest)
sosUOM(sosResult(obs.temp.latest))

################################################################################
# adding filename attribute testing
weathersos <- SOS(SosExampleServices()[[1]])

# o_4995049, o_4995048
obs <- getObservationById(weathersos, "o_4995049", verbose = TRUE,
		saveOriginal = TRUE)
attributes(obs)
attributes(obs)[[sosAttributeFileName]]
# works!

obs2 <- getObservation(weathersos, offering = sosOfferings(weathersos)[[1]],
#		verbose = TRUE,
		latest = TRUE,
		saveOriginal = "testObservation"
)
attributes(obs2)
# works!

sensor <- describeSensor(weathersos, sosProcedures(weathersos)[[1]][[1]],
		saveOriginal = "mySensor")
attributes(sensor)
# works!


################################################################################
# Google CRS
require("rgdal", quietly = TRUE)
CRS("+init=epsg:4979")
# raises error

?try

lala <- function() {
	.temp <- NULL
#	try(.temp <- CRS("+init=epsg:4979"))
#	print("done 1")
#	print(.temp)
	
	tryCatch({
				.temp <- CRS("+init=epsg:4979")
			}, error = function(ex) {
				cat("An error was detected: ", toString(ex))
			}, finally = {
				cat("Releasing resources...");
				cat("done\n");
			})
	print("done 1")
	print(.temp)
}
lala()

sosGetCRS("epsg:1234")


################################################################################
# Ben tips about 3d plots in R
library(gstat)
library(rgl)
library(mgcv)
library(car)

data(meuse)

scatter3d(meuse$dist, meuse$lead, meuse$elev, fit="linear", 
		residuals=TRUE, groups=meuse$ffreq, parallel=TRUE, bg="white", 
		axis.scales=TRUE, grid=TRUE, ellipsoid=FALSE, xlab="dist", ylab="lead", 
		zlab="elev")


################################################################################
# Try StatET 2.0's new visual debugger 
testfunction <- function(x, y) {
	for (i in seq(from = 1, to = x)) {
		print(i + y)
		
		if(i %% 2 == 1) { # odd number
			var <- othertestfunction(i, y)
			
			print(var)
		}
	}
}

othertestfunction <- function(x, z) {
	return(x - z)
}

# loading sos4R via library does not work!

################################################################################
# test parsing multiline file for shorte lines in example section of .Rd files
library("XML")
er.doc <- xmlParseDoc(paste0("<ows:ExceptionReport ", 
	"xmlns:ows=\"http://www.opengis.net/ows/1.1\" ",
	"xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" ",
	"version=\"1.0.0\" ",
	"xsi:schemaLocation=",
	"\"http://schemas.opengis.net/ows/1.1.0/owsExceptionReport.xsd\">",
	"<ows:Exception exceptionCode=\"VersionNegotiationFailed\" ",
	"locator=\"AcceptVersions\">",
	"<ows:ExceptionText>The parameter 'AcceptVersions' does not contain the ",
	"version of this SOS: '1.0.0'</ows:ExceptionText>",
	"</ows:Exception>",
	"</ows:ExceptionReport>"))
er.doc
