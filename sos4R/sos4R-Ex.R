pkgname <- "sos4R"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
options(pager = "console")
library('sos4R')

base::assign(".oldSearch", base::search(), pos = 'CheckExEnv')
cleanEx()
nameEx("Constants")
### * Constants

flush(stderr()); flush(stdout())

### Name: Constants
### Title: Constants in sos4R
### Aliases: Constants gmlBeginName gmlBeginPositionName gmlBoundedByName
###   gmlDescriptionName gmlDurationName gmlEndName gmlEndPositionName
###   gmlEnvelopeName gmlFeatureCollectionName gmlFeatureMemberName
###   gmlLowerCornerName gmlNameName gmlNamespacePrefix gmlPointName
###   gmlPosName gmlRelatedTimeName gmlTimeInstantName gmlTimeIntervalName
###   gmlTimeLengthName gmlTimePeriodName gmlTimePositionName
###   gmlUpperCornerName ogcBBOXName ogcArithmeticOperatorsName
###   ogcComparisonOpBetweenName ogcComparisonOpEqualToName
###   ogcComparisonOpGreaterThanName
###   ogcComparisonOpGreaterThanOrEqualToName ogcComparisonOpLessThenName
###   ogcComparisonOpLessThanOrEqualToName ogcComparisonOpIsLikeName
###   ogcComparisonOpIsNotEqualTo ogcComparisonOpIsNull ogcContainsName
###   ogcEIDName ogcFIDName ogcGeometryOperandName ogcGeometryOperandsName
###   ogcIdCapabilities ogcGeometryOperandEnvelopeName
###   ogcGeometryOperandPolygonName ogcGeometryOperandPointName
###   ogcGeometryOperandLineStringName ogcIntersectsName ogcLiteralName
###   ogcLogicalOperatorsName ogcNamespacePrefix ogcOverlapsName
###   ogcPropertyNameName ogcScalarCapabilitiesName
###   ogcSpatialCapabilitiesName ogcSpatialOperatorName
###   ogcSpatialOperatorsName ogcSpatialOpBBOXName ogcSpatialOpContainsName
###   ogcSpatialOpIntersectsName ogcSpatialOpOverlapsName
###   ogcSpatialOpBeyondName ogcSpatialOpCrossesName
###   ogcSpatialOpDWithinName ogcSpatialOpDisjointName
###   ogcSpatialOpEqualsName ogcSpatialOpTouchesName ogcSpatialOpWithinName
###   ogcTempOpTMAfterName ogcTempOpTMBeforeName ogcTempOpTMBeginsName
###   ogcTempOpTMBegunByName ogcTempOpTMContainsName ogcTempOpTMDuringName
###   ogcTempOpTMEndedByName ogcTempOpTMEndsName ogcTempOpTMEqualsName
###   ogcTempOpTMMeetsName ogcTempOpTMMetByName ogcTempOpTMOverlapsName
###   ogcTempOpTTMOverlappedBy ogcTemporalCapabilitiesName
###   ogcTemporalOperandsName ogcTemporalOperandName
###   ogcTemporalOperatorsName ogcTemporalOperatorName omMeasurementName
###   omMemberName omObservationName omObservationCollectionName
###   omFeatureOfInterestName omProcedureName omObservedPropertyName
###   omResultTimeName omSamplingTimeName omResultName
###   omCategoryObservationName omCountObservationName
###   omTruthObservationName omGeometryObservationName
###   omTemporalObservationName omComplexObservationName
###   saSamplingPointName saSamplingSurface saPositionName
###   saSampledFeatureName saSamplingTimeName owsServiceIdentificationName
###   owsTitleName owsAbstractName owsKeywordsName owsKeywordName
###   owsServiceTypeName owsServiceTypeVersionName owsFeesName
###   owsAccessConstraintsName owsServiceProviderName
###   owsOperationsMetadataName owsOperationName owsDCPName owsHTTPName
###   owsGetName owsPostName owsParameterName owsAllowedValuesName
###   owsValueName owsAnyValueName owsRangeName owsMinimumValueName
###   owsMaximumValueName owsSpacingName owsConstraintName owsMetadataName
###   owsExceptionName owsExceptionTextName owsProfileName
###   owsProviderNameName owsProviderSiteName owsServiceContactName
###   saSamplingPointName saSamplingSurface saPositionName
###   saSampledFeatureName saSamplingTimeName sosService sosNamespacePrefix
###   sosGetCapabilitiesName sosDescribeSensorName sosGetObservationName
###   sosGetObservationByIdName owsExceptionReportName
###   sosGetFeatureOfInterestName sosIntendedApplicationName sosTimeName
###   sosProcedureName sosObservedPropertyName sosFeatureOfInterestName
###   sosResultModelName sosResponseFormatName sosResponseModeName
###   sosObservationOfferingName sosObservationOfferingListName
###   sosContentsName sosFilterCapabilitiesName sosCapabilitiesName
###   sosEventTimeName sosEventTimeLatestValue sosObjectIDName
###   sosResultName sweCompositePhenomenonName sweBaseName sweComponentName
###   sweDataArrayName sweElementTypeName sweSimpleDataRecordName
###   sweDataRecordName sweFieldName sweTimeName sweQuantityName
###   sweCategoryName sweBooleanName sweCountName sweEncodingName
###   sweTextBlockName sweValuesName sweValueName sweCodeSpaceName
###   sweTextName sweUomName sweCoordinateName sweLocationName
###   swePositionName sweVectorName xmlInternalDocumentName xmlTextNodeName
###   OwsExceptionsData ogcComparisonOperatorsName ogcTempOpTMOverlappedBy
###   owsNamespacePrefix sosDescribeFeatureTypeName
###   sosDescribeObservationTypeName sosDescribeResultModelName
###   sosGetFeatureOfInterestTimeName sosGetResultName
###   sosInsertObservationName sosRegisterSensorName mimeTypeCSV mimeTypeOM
###   mimeTypeSML mimeTypeXML mimeSubtypeOM smlSensorMLName
###   sosAttributeFileName
### Keywords: constants XML

### ** Examples


# example constants
sosNamespacePrefix
gmlNameName
sweUomName

# Data frame holding OWS exception code information
OwsExceptionsData()



cleanEx()
nameEx("Defaults")
### * Defaults

flush(stderr()); flush(stdout())

### Name: Defaults
### Title: Default Parameter Settings and Handling Functions
### Aliases: sosDefault Defaults sosDefaultCharacterEncoding
###   SosDefaultConnectionMethod sosDefaultDescribeSensorOutputFormat
###   sosDefaultGetCapAcceptFormats sosDefaultGetCapOwsVersion
###   sosDefaultGetCapSections sosDefaultGetObsResponseFormat
###   sosDefaultSpatialOpPropertyName sosDefaultTempOpPropertyName
###   sosDefaultTemporalOperator sosDefaultTimeFormat
###   sosDefaultFilenameTimeFormat sosDefaultColumnNameFeatureIdentifier
###   sosDefaultColumnNameLat sosDefaultColumnNameLon
###   sosDefaultColumnNameSRS sosDefaultColorPalette
###   sosDefaultReferenceFrameSensorDescription SosParsingFunctions
###   SosEncodingFunctions SosDisabledParsers
###   SosDataFieldConvertingFunctions SosExampleServices SosDefaults
###   SosResetParsingFunctions
### Keywords: misc

### ** Examples

# simple default values
show(sosDefaultCharacterEncoding)
show(sosDefaultDescribeSensorOutputFormat)
show(sosDefaultGetCapAcceptFormats)
show(sosDefaultGetCapOwsVersion)
show(sosDefaultGetCapSections)
show(sosDefaultGetObsResponseFormat)
show(sosDefaultSpatialOpPropertyName)
show(sosDefaultTempOpPropertyName)
show(sosDefaultTemporalOperator)
show(sosDefaultTimeFormat)
SosDefaultConnectionMethod()

## Not run: 
##D # usage of defaults in construction method for SOS class
##D sos <- SOS("http://mysos.com/sos", method = SosDefaultConnectionMethod(),
##D 		timeFormat = sosDefaultTimeFormat)
## End(Not run)

# functions to disable all parsing
SosDisabledParsers()

# Replace a parsing function
myER <- function(xml) {
	return("EXCEPTION!!!11")
}
SosParsingFunctions("ExceptionReport" = myER)

# use inclusion and exclusion, important: even the just added function needs to be included manually!
SosParsingFunctions("ExceptionReport" = myER,
	include = c("GetObservation", "DescribeSensor", "ExceptionReport"))
SosParsingFunctions(exclude = c("GetObservation", "DescribeSensor"))

## Not run: 
##D # Replace an encoding function
##D myEncoding <- function(object, v) {
##D 	return(str(object))
##D }
##D 
##D sos = SOS(url = "http://mysos.com/sos",
##D 		encoders = SosEncodingFunctions("POST" = myPostEncoding))
##D 
##D # Use custom converting function and connection method. This mechanism works the same for encoders and decoders.
##D myConverters <- SosDataFieldConvertingFunctions(
##D 	"myNumericUnit" = sosConvertDouble,
##D mySos <- SOS(sos.url, method = "GET", dataFieldConverters = myConverters)
##D sosDataFieldConverters(mySos)
##D 
##D # inspecting XML using dummy parsing function
##D sos = SOS(url = "http://mysos.com/sos", parsers = SosDisabledParsers)
##D describeSensor(sos, sosProcedures(sos)[[1]])
## End(Not run)

# a list of example services
SosExampleServices()

# a named list of all defaults
SosDefaults()

# replace the parsing functions with the default ones
## Not run: 
##D sos <- SosResetParsingFunctions(sos)
## End(Not run)




cleanEx()
nameEx("DescribeSensor")
### * DescribeSensor

flush(stderr()); flush(stdout())

### Name: DescribeSensor
### Title: Class and Construction Function for "SosDescribeSensor"
### Aliases: DescribeSensor SosDescribeSensor SosDescribeSensor-class
###   show,SosDescribeSensor-method print,SosDescribeSensor-method
###   toString,SosDescribeSensor-method
### Keywords: classes

### ** Examples

showClass("SosDescribeSensor")

# example for construction function
describeSensorRequest <- SosDescribeSensor(service = "SOS", version = "1.0.0",
	procedure = "urn:procedure:42", outputFormat = "text/xml")
print(describeSensorRequest)

# encode the request in XML
encodeRequestXML(describeSensorRequest)




cleanEx()
nameEx("GML")
### * GML

flush(stderr()); flush(stdout())

### Name: GmlDirectPosition-class
### Title: Classes and Construction Functions from the GML Namespace
### Aliases: GmlDirectPosition-class GmlDirectPositionOrNULL-class
###   GmlEnvelope-class GmlFeature-class GmlFeatureCollection-class
###   GmlFeatureOrNULL-class GmlFeatureProperty-class
###   GmlFeatureOrGmlFeaturePropertyOrNULL-class GmlGeometry-class
###   GmlLineString-class GmlPoint-class GmlPointOrNULL-class
###   GmlPointProperty-class GmlPolygon-class
###   GmlTimeGeometricPrimitive-class GmlTimeInstant-class
###   GmlTimeInstantOrNULL-class GmlTimeInstantProperty-class
###   GmlTimeInstantPropertyOrNULL-class GmlTimeInterval-class
###   GmlTimeIntervalOrNULL-class GmlTimeObject-class
###   GmlTimeObjectOrNULL-class GmlTimePeriod-class GmlTimePosition-class
###   GmlTimePositionOrNULL-class GmlTimePrimitive-class
###   show,GmlDirectPosition-method show,GmlEnvelope-method
###   show,GmlFeatureCollection-method show,GmlFeatureProperty-method
###   show,GmlGeometry-method show,GmlPoint-method
###   show,GmlPointProperty-method show,GmlTimeInstant-method
###   show,GmlTimeInterval-method show,GmlTimePeriod-method
###   show,GmlTimePosition-method print,GmlDirectPosition-method
###   print,GmlEnvelope-method print,GmlFeatureCollection-method
###   print,GmlFeatureProperty-method print,GmlGeometry-method
###   print,GmlPoint-method print,GmlPointProperty-method
###   print,GmlTimeInstant-method print,GmlTimeInstantProperty-method
###   print,GmlTimeInterval-method print,GmlTimePeriod-method
###   print,GmlTimePosition-method toString,GmlDirectPosition-method
###   toString,GmlEnvelope-method toString,GmlFeatureCollection-method
###   toString,GmlFeatureProperty-method toString,GmlGeometry-method
###   toString,GmlPoint-method toString,GmlPointProperty-method
###   toString,GmlTimeInstant-method toString,GmlTimeInstantProperty-method
###   toString,GmlTimeInterval-method toString,GmlTimePeriod-method
###   toString,GmlTimePosition-method GmlDirectPosition
###   GmlDirectPositionLatLon GmlEnvelope GmlFeatureCollection GmlPoint
###   GmlPointProperty GmlFeatureProperty GmlTimeInstant
###   GmlTimeInstantProperty GmlTimeInterval GmlTimePeriod GmlTimePosition
###   GmlMeasure-class print,GmlMeasure-method toString,GmlMeasure-method
###   GmlMeasure show,GmlMeasure-method
###   sosCoordinates,GmlDirectPosition-method
###   sosCoordinates,GmlFeatureCollection-method
###   sosCoordinates,GmlFeatureProperty-method
###   sosCoordinates,GmlPoint-method sosCoordinates,GmlPointProperty-method
###   sosId,GmlFeature-method sosSrsName,GmlDirectPosition-method
###   sosSrsName,GmlPoint-method sosFeatureIds,GmlFeatureCollection-method
###   sosFeatureIds,GmlFeatureProperty-method
###   sosFeaturesOfInterest,GmlFeatureCollection-method
###   sosTime,GmlTimeInstant-method sosTime,GmlTimeInstantProperty-method
###   sosTime,GmlTimePeriod-method sosTime,GmlTimePosition-method
###   summary.GmlTimePeriod print.summary.GmlTimePeriod
###   sosUOM,GmlMeasure-method
### Keywords: classes utilities

### ** Examples

showClass("GmlDirectPosition")
showClass("GmlEnvelope")
showClass("GmlFeature")
showClass("GmlFeatureCollection")
showClass("GmlFeatureOrNULL")
showClass("GmlFeatureProperty")
showClass("GmlGeometry")
showClass("GmlLineString")
showClass("GmlPoint")
showClass("GmlPointProperty")
showClass("GmlPolygon")
showClass("GmlTimeGeometricPrimitive")
showClass("GmlTimeInstant")
showClass("GmlTimeInstantOrNULL")
showClass("GmlTimeInstantProperty")
showClass("GmlTimeInstantPropertyOrNULL")
showClass("GmlTimeInterval")
showClass("GmlTimeIntervalOrNULL")
showClass("GmlTimeObject")
showClass("GmlTimeObjectOrNULL")
showClass("GmlTimePeriod")
showClass("GmlTimePosition")
showClass("GmlTimePositionOrNULL")
showClass("GmlTimePrimitive")

# create direct position
pos1 <- GmlDirectPosition(pos = "7.0 52.0")
show(pos1)

# create envelope
env1 <- GmlEnvelope(upperCorner = pos1, lowerCorner = GmlDirectPosition("6.0 51.0"))
print(env1)

# wrap elements in feature collection
GmlFeatureCollection(id = "001", featureMembers=list(pos1, env1))

# create point with ID
point1 <- GmlPoint(pos = pos1, id = "002")

# create point properties
GmlPointProperty(href = "http://link.to/point")
GmlPointProperty(point = point1)

# time interval of one day
GmlTimeInterval(interval = "1", unit = "d")

# referenced feature
GmlFeatureProperty(href = "http://link.to/feature")

# create a time position and wrap it into a time instant
timePos1 <- GmlTimePosition(time = as.POSIXct("2010-01-01"))

# create direct or referenced time instant
timeInst1 <- GmlTimeInstant(timePosition = timePos1)
timeInst1

GmlTimeInstantProperty(href = "http://link.to/timeInstant")

# create different variants of time periods
# one hour with time positions
GmlTimePeriod(beginPosition = timePos1, endPosition = GmlTimePosition(time = timePos1@time+3600))

# one week backwards from now
aWeekAgo <- GmlTimeInstantProperty(time = GmlTimeInstant(time = GmlTimePosition(time = Sys.time()-(3600*24*7))))
now <- GmlTimeInstantProperty(time = GmlTimeInstant(time = GmlTimePosition(time = Sys.time())))
GmlTimePeriod(begin = aWeekAgo, end = now)




cleanEx()
nameEx("GetObservation")
### * GetObservation

flush(stderr()); flush(stdout())

### Name: GetObservation
### Title: GetObservation and GetObservationById Request Objects
### Aliases: GetObservation SosGetObservation SosGetObservation-class
###   show,SosGetObservation-method toString,SosGetObservation-method
###   print,SosGetObservation-method GetObservationById
###   SosGetObservationById SosGetObservationById-class
###   show,SosGetObservationById-method print,SosGetObservationById-method
###   toString,SosGetObservationById-method
### Keywords: classes utitlities

### ** Examples

showClass("SosGetObservation")
showClass("SosGetObservationById")

observationRequest <- SosGetObservation(service = "SOS", version = "1.0.0", offering = "temperatures", observedProperty = list("urn:property:AirTemperature"), responseFormat = "text/xml;subtype=&quot;om/1.0.0&quot;")
print(observationRequest)

observationByIdRequest <- SosGetObservationById(service = "SOS", version = "1.0.0", observationId = "o_12345", responseFormat = "text/xml;subtype=&quot;om/1.0.0&quot;")
print(observationByIdRequest)

## Not run: 
##D sos <- SOS("http://mysos.net/sos")
##D encodeXML(observationByIdRequest, sos = sos)
## End(Not run)




cleanEx()
nameEx("KML")
### * KML

flush(stderr()); flush(stdout())

### Name: KML
### Title: Methods for the Namespace kml
### Aliases: parseKML mimeTypeKML kmlName kml
### Keywords: methods misc

### ** Examples

#



cleanEx()
nameEx("OGC")
### * OGC

flush(stderr()); flush(stdout())

### Name: OGC
### Title: Classes and Construction Functions for the OGC Namespace
### Aliases: OGC ogc OgcBBOX-class show,OgcBBOX-method OgcBinarySpatialOp
###   OgcBinarySpatialOp-class OgcBinaryTemporalOp
###   OgcBinaryTemporalOp-class OgcBinaryTemporalOpOrNULL-class
###   OgcComparisonOps OgcComparisonOps-class OgcContains-class
###   show,OgcContains-method OgcIntersects-class show,OgcIntersects-method
###   OgcOverlaps-class show,OgcOverlaps-method OgcSpatialOps
###   OgcSpatialOps-class OgcSpatialOpsOrNULL-class OgcBBOX OgcContains
###   OgcIntersects OgcOverlaps print,OgcBBOX-method
###   print,OgcContains-method print,OgcIntersects-method
###   print,OgcOverlaps-method toString,OgcBBOX-method
###   toString,OgcContains-method toString,OgcIntersects-method
###   toString,OgcOverlaps-method
### Keywords: classes utilities

### ** Examples

showClass("OgcBBOX")
showClass("OgcBinarySpatialOp")
showClass("OgcBinaryTemporalOp")
showClass("OgcBinaryTemporalOpOrNULL")
showClass("OgcComparisonOps")
showClass("OgcContains")
showClass("OgcOverlaps")
showClass("OgcSpatialOps")
showClass("OgcSpatialOpsOrNULL")

# TBD examples for construction functions



cleanEx()
nameEx("OWS")
### * OWS

flush(stderr()); flush(stdout())

### Name: OWS
### Title: Classes and Construction Functions for Elements of the OWS
###   Namespace
### Aliases: OwsCapabilities_1.1.0-class OwsCapabilities_2.0.0-class
###   OwsCapabilities-class OwsContents-class OwsContentsOrNULL-class
###   OwsException-class OwsExceptionReport OwsExceptionReport-class
###   OwsGetCapabilities_1.1.0-class OwsGetCapabilities_2.0.0-class
###   OwsGetCapabilities-class OwsOperation-class
###   OwsOperationsMetadata-class OwsOperationsMetadataOrNULL-class
###   OwsRange-class OwsServiceIdentification-class
###   OwsServiceIdentificationOrNULL-class OwsServiceOperation-class
###   OwsServiceProvider-class OwsServiceProviderOrNULL-class
###   print,OwsCapabilities-method print,OwsCapabilities_1.1.0-method
###   print,OwsCapabilities_2.0.0-method print,OwsContents-method
###   print,OwsException-method print,OwsExceptionReport-method
###   print,OwsGetCapabilities-method print,OwsGetCapabilities_1.1.0-method
###   print,OwsGetCapabilities_2.0.0-method print,OwsOperation-method
###   print,OwsOperationsMetadata-method print,OwsRange-method
###   print,OwsServiceIdentification-method
###   print,OwsServiceOperation-method print,OwsServiceProvider-method
###   print.summary.OwsRange show,OwsCapabilities_1.1.0-method
###   show,OwsCapabilities_2.0.0-method show,OwsCapabilities-method
###   show,OwsContents-method show,OwsException-method
###   show,OwsExceptionReport-method show,OwsGetCapabilities_1.1.0-method
###   show,OwsGetCapabilities_2.0.0-method show,OwsGetCapabilities-method
###   show,OwsOperation-method show,OwsOperationsMetadata-method
###   show,OwsRange-method show,OwsServiceIdentification-method
###   show,OwsServiceOperation-method show,OwsServiceProvider-method
###   toString,OwsCapabilities-method toString,OwsCapabilities_1.1.0-method
###   toString,OwsCapabilities_2.0.0-method toString,OwsContents-method
###   toString,OwsException-method toString,OwsExceptionReport-method
###   toString,OwsGetCapabilities-method
###   toString,OwsGetCapabilities_1.1.0-method
###   toString,OwsGetCapabilities_2.0.0-method toString,OwsOperation-method
###   toString,OwsOperationsMetadata-method toString,OwsRange-method
###   toString,OwsServiceIdentification-method
###   toString,OwsServiceOperation-method
###   toString,OwsServiceProvider-method OwsCapabilities OwsContents
###   OwsException OwsGetCapabilities OwsOperation OwsOperationsMetadata
###   OwsRange OwsServiceIdentification OwsServiceProvider
###   sosResult,OwsExceptionReport-method
###   sosTitle,OwsServiceIdentification-method
###   sosAbstract,OwsServiceIdentification-method summary.OwsRange
### Keywords: classes utilities

### ** Examples

showClass("OwsCapabilities_1.1.0")
showClass("OwsCapabilities_2.0.0")
showClass("OwsCapabilities")
showClass("OwsContents")
showClass("OwsContentsOrNULL")
showClass("OwsException")
showClass("OwsExceptionReport")
showClass("OwsGetCapabilities_1.1.0")
showClass("OwsGetCapabilities_2.0.0")
showClass("OwsGetCapabilities")
showClass("OwsOperation")
showClass("OwsOperationsMetadata")
showClass("OwsRange")
showClass("OwsServiceIdentification")
showClass("OwsServiceIdentificationOrNULL")
showClass("OwsServiceOperation")
showClass("OwsServiceProvider")
showClass("OwsServiceProviderOrNULL")

# TBD examples for construction functions




cleanEx()
nameEx("OmMeasurement")
### * OmMeasurement

flush(stderr()); flush(stdout())

### Name: OmMeasurement
### Title: Class and Construction Function for om:Measurement Elements
### Aliases: OmMeasurement OmMeasurement-class show,OmMeasurement-method
###   print,OmMeasurement-method toString,OmMeasurement-method
###   sosResult,OmMeasurement-method sosProcedures,OmMeasurement-method
###   sosFeatureIds,OmMeasurement-method
###   sosFeaturesOfInterest,OmMeasurement-method
###   sosGetCRS,OmMeasurement-method names.OmMeasurement
###   as.data.frame.OmMeasurement as.SpatialPointsDataFrame.OmMeasurement
###   sosUOM,OmMeasurement-method
### Keywords: classes

### ** Examples

showClass("OmMeasurement")

# TBD examples for construction function

# TBD examples for sosResult



cleanEx()
nameEx("OmObservation")
### * OmObservation

flush(stderr()); flush(stdout())

### Name: OmObservation-class
### Title: Classes for om:Observation Elements
### Aliases: OmObservation OmObservation-class OmObservationOrNULL-class
###   show,OmObservation-method print,OmObservation-method
###   print,OmObservationProperty-method toString,OmObservation-method
###   toString,OmObservationProperty-method OmObservationProperty-class
###   show,OmObservationProperty-method OmObservationProperty
###   sosProcedures,OmObservation-method names.OmObservation sosResult
###   sosResult,list-method sosResult,OmObservationProperty-method
###   sosResult,OmObservation-method as.data.frame.OmObservation
###   sosCoordinates,OmObservation-method
###   sosFeaturesOfInterest,OmObservation-method
###   sosFeatureIds,OmObservation-method
###   sosObservedProperties,OmObservation-method
###   sosGetCRS,OmObservation-method
###   as.SpatialPointsDataFrame.OmObservation sosUOM,OmObservation-method
###   print.summary.OmObservation summary.OmObservation
### Keywords: classes

### ** Examples

showClass("OmObservation")
showClass("OmObservationProperty")
showClass("OmObservationOrNULL")

# TBD examples for construction methods
OmObservationProperty(href = "http://link.to/myObservation")

# get result from an observation
## Not run: 
##D result <- observation@result
##D 
##D # the accessor method also works with lists of observations
##D result <- sosResult(observation)
##D resultList <- sosResult(observationList)
## End(Not run)




cleanEx()
nameEx("OmObservationCollection")
### * OmObservationCollection

flush(stderr()); flush(stdout())

### Name: OmObservationCollection
### Title: Class "OmObservationCollection"
### Aliases: OmObservationCollection OmObservationCollection-class
###   length,OmObservationCollection-method
###   show,OmObservationCollection-method
###   sosResult,OmObservationCollection-method
###   print,OmObservationCollection-method
###   toString,OmObservationCollection-method
###   [,OmObservationCollection-method
###   [[,OmObservationCollection,ANY,missing-method
###   as.list.OmObservationCollection length.OmObservationCollection
###   names.OmObservationCollection
###   sosBoundedBy,OmObservationCollection-method
###   sosCoordinates,OmObservationCollection-method
###   sosProcedures,OmObservationCollection-method
###   sosFeatureIds,OmObservationCollection-method
###   sosObservedProperties,OmObservationCollection
###   sosObservedProperties,OmObservationCollection-method
###   sosFeaturesOfInterest,OmObservationCollection
###   sosFeaturesOfInterest,OmObservationCollection-method
###   as.SpatialPointsDataFrame.OmObservationCollection
###   sosGetCRS,OmObservationCollection-method
###   sosUOM,OmObservationCollection-method
###   print.summary.OmObservationCollection summary.OmObservationCollection
### Keywords: classes

### ** Examples

showClass("OmObservationCollection")




cleanEx()
nameEx("SA")
### * SA

flush(stderr()); flush(stdout())

### Name: SA
### Title: Classes of the Namespace sa
### Aliases: sa 'sampling features' SaSamplingPoint SaSamplingPoint-class
###   show,SaSamplingPoint-method SaSamplingSurface SaSamplingSurface-class
###   show,SaSamplingSurface-method print,SaSamplingPoint-method
###   print,SaSamplingSurface-method toString,SaSamplingPoint-method
###   toString,SaSamplingSurface-method SaSamplingPoint
###   sosCoordinates,SaSamplingPoint-method
###   sosFeatureIds,SaSamplingPoint-method
### Keywords: classes

### ** Examples

showClass("SaSamplingPoint")

# create sampling point
SaSamplingPoint(sampledFeatures = list("feature1", "feature2"), position = GmlPointProperty(href = "http://link.to/point"))




cleanEx()
nameEx("SML")
### * SML

flush(stderr()); flush(stdout())

### Name: SML
### Title: Classes of the Namespace sml
### Aliases: SensorML-class show,SensorML-method print,SensorML-method
###   toString,SensorML-method SensorML sml sosId,SensorML-method
###   sosName,SensorML-method sosAbstract,SensorML-method
###   sosCoordinates,SensorML-method sosBoundedBy,SensorML-method
###   sosGetCRS,SensorML-method as.SensorML.SpatialPointsDataFrame
###   plot.SensorML plot,SensorML,missing-method
### Keywords: classes

### ** Examples

showClass("SensorML")

## Not run: 
##D weathersos <- SOS("http://v-swe.uni-muenster.de:8080/WeatherSOS/sos")
##D proc1 <- sosProcedures(weathersos)[[1]][[1]]
##D proc1.descr <- describeSensor(weathersos, proc1, verbose = TRUE)
##D plot(proc1.descr)
##D class(proc1.descr)
##D print(proc1.descr
## End(Not run)




cleanEx()
nameEx("SOS")
### * SOS

flush(stderr()); flush(stdout())

### Name: SOS
### Title: Class, and Construction and Accessor Functions for "SOS"
### Aliases: SOS SOS-class show,SOS-method print,SOS-method
###   toString,SOS-method SOS_1.0.0 SOS_1.0.0-class show,SOS_1.0.0-method
###   print,SOS_1.0.0-method toString,SOS_1.0.0-method
###   sosCapabilitiesDocumentOriginal
###   sosCapabilitiesDocumentOriginal,SOS-method sosCaps sosCaps-methods
###   sosCaps,SOS-method sosContents sosContents-methods
###   sosContents,SOS-method sosDataFieldConverters
###   sosDataFieldConverters-methods sosDataFieldConverters,SOS-method
###   sosTime sosTime-methods sosTime,SOS-method sosTime,list-method
###   sosFilter_Capabilities sosFilter_Capabilities-methods
###   sosFilter_Capabilities,SOS-method sosFeaturesOfInterest
###   sosFeaturesOfInterest-methods sosFeaturesOfInterest,SOS-method
###   sosFeaturesOfInterest,SOS,character-method
###   sosFeaturesOfInterest,SosObservationOffering-method sosMethod
###   sosMethod-methods sosMethod,SOS-method sosMethod,SOS_1.0.0-method
###   sosObservedProperties sosObservedProperties-methods
###   sosObservedProperties,SOS-method
###   sosObservedProperties,SosObservationOffering-method sosOfferingIds
###   sosOfferingIds-methods sosOfferingIds,SOS-method sosOfferings
###   sosOfferings-methods sosOfferings,SOS-method
###   sosOfferings,SOS,character-method sosOperation sosOperation-methods
###   sosOperation,SOS,character-method sosOperationsMetadata
###   sosOperationsMetadata-methods sosOperationsMetadata,SOS-method
###   sosParsers sosParsers-methods sosParsers,SOS-method sosProcedures
###   sosProcedures-methods sosProcedures,SOS-method
###   sosProcedures,list-method sosProcedures,SosObservationOffering-method
###   sosResponseFormats sosResponseFormats-methods
###   sosResponseFormats,SOS-method sosResponseFormats,OwsOperation-method
###   sosResponseFormats,SosObservationOffering-method sosResponseMode
###   sosResponseMode-methods sosResponseMode,SOS-method
###   sosResponseMode,OwsOperation-method
###   sosResponseMode,SosObservationOffering-method sosResultModels
###   sosResultModels-methods sosResultModels,SOS-method
###   sosResultModels,OwsOperation-method sosResult,character-method
###   sosServiceIdentification sosServiceIdentification-methods
###   sosServiceIdentification,SOS-method sosServiceProvider
###   sosServiceProvider-methods sosServiceProvider,SOS-method sosSrsName
###   sosSrsName-methods sosSrsName,SOS-method sosTimeFormat
###   sosTimeFormat-methods sosTimeFormat,SOS-method sosUrl sosUrl-methods
###   sosUrl,SOS-method sosUrl,SOS_1.0.0-method sosVersion
###   sosVersion-methods sosVersion,SOS-method sosExceptionCodeMeaning
###   sosExceptionCodeMeaning,character-method sosBoundedBy
###   sosBoundedBy-method sosBoundedBy,list-method sosCoordinates
###   sosCoordinates-method sosCoordinates,list-method
###   sosCoordinates,SosObservationOffering-method sosId sosId-method
###   sosId,list-method sosSrsName sosSrsName-method sosFeatureIds
###   sosFeatureIds-method sosFeatureIds,list-method
###   sosObservedProperties,list-method sosFeaturesOfInterest,list-method
###   sosGetCRS sosGetCRS-method sosGetCRS,SOS-method sosGetCRS,list-method
###   sosGetCRS,character-method sosGetCRS,SosObservationOffering-method
###   sosName sosName-method sosName,list-method
###   sosName,OwsServiceProvider-method sosName,OwsGetCapabilities-method
###   sosName,OwsOperation-method sosName,SosDescribeSensor-method
###   sosName,SosGetObservation-method sosName,SosGetObservationById-method
###   sosTitle sosTitle-method sosTitle,SOS-method sosAbstract
###   sosAbstract-method sosAbstract,SOS-method sosEncoders
###   sosEncoders,SOS-method sosOperations sosOperations,SOS-method
###   sosOperations,OwsCapabilities-method
###   sosOperations,SosCapabilities_1.0.0-method sosGetDCP
###   sosGetDCP,SOS,character-method sosSwitchCoordinates
###   sosSwitchCoordinates-method sosSwitchCoordinates,SOS-method plot.SOS
###   plot,SOS,missing-method plot.SosObservationOffering
###   plot,SosObservationOffering,missing-method
###   as.SosObservationOffering.SpatialPolygons print.summary.SOS
###   print.summary.SosObservationOffering summary.SosObservationOffering
###   summary.SOS sosResult,data.frame-method sosUOM sosUOM,list-method
###   sosUOM,data.frame-method sosCapabilitiesUrl sosCapabilitiesUrl-method
###   sosCapabilitiesUrl,SOS-method
### Keywords: classes

### ** Examples

showClass("SOS")

## Not run: 
##D # create a SOS connection
##D mysos <- SOS(url = "http://mysos.org/sos")
##D 
##D # create a SOS connetion with a specific connection method and time format
##D mysos <- SOS(url = "http://mysos.org/sos",
##D     method = "GET", timeFormat = "##D 
##D 
##D # turn on verbose output for all methods and functions
##D SOS(url = "http://mysos.org/sos", verboseOutput = TRUE)
##D 
##D # get the meaning of an exception code
##D sosExceptionCodeMeaning(ex@exceptionCode)
##D 
##D # create a CRS object from a URN CRS string
##D sosGetCRS("urn:ogc:def:crs:EPSG:4217")
##D 
##D # create the URL to a GET request for GetCapabilities
##D sosCapabilitiesUrl(mysos)
## End(Not run)






cleanEx()
nameEx("SWE")
### * SWE

flush(stderr()); flush(stdout())

### Name: SWE
### Title: Classes and Construction Functions for the SWE Namespace
### Aliases: SweCompositePhenomenon-class
###   show,SweCompositePhenomenon-method SwePhenomenon-class
###   show,SwePhenomenon-method SwePhenomenonOrNULL-class
###   SwePhenomenonProperty-class show,SwePhenomenonProperty-method
###   SwePhenomenonPropertyOrNULL-class SweTextBlock-class
###   show,SweTextBlock-method print,SweCompositePhenomenon-method
###   print,SwePhenomenon-method print,SwePhenomenonProperty-method
###   print,SweTextBlock-method toString,SweCompositePhenomenon-method
###   toString,SwePhenomenon-method toString,SwePhenomenonProperty-method
###   toString,SweTextBlock-method SweCompositePhenomenon SwePhenomenon
###   SwePhenomenonProperty SweTextBlock
###   sosObservedProperties,SweCompositePhenomenon-method
###   sosObservedProperties,SwePhenomenonProperty-method
### Keywords: classes

### ** Examples

showClass("SweCompositePhenomenon")
showClass("SwePhenomenon")
showClass("SwePhenomenonProperty")
showClass("SwePhenomenonPropertyOrNULL")
showClass("SweTextBlock")


# TBD examples for construction functions




cleanEx()
nameEx("SosBindings")
### * SosBindings

flush(stderr()); flush(stdout())

### Name: SosBindings
### Title: Bindings and Connecition Methods of OGC Sensor Observation
###   Service
### Aliases: 'connection methods' bindings SosBinding SosBindings GET POST
###   SOAP HTTP
### Keywords: constants XML

### ** Examples


# HTTP connection methods supported by this sos4R implementation
supported <- SosSupportedConnectionMethods()
supported

## Not run: 
##D sos <- SOS("http://sosurl.org/", method = "GET")
## End(Not run)




cleanEx()
nameEx("SosCapabilities")
### * SosCapabilities

flush(stderr()); flush(stdout())

### Name: SosCapabilities_1.0.0-class
### Title: Class and Construction Function for "SosCapabilities_1.0.0"
### Aliases: SosCapabilities_1.0.0-class SosCapabilities
###   SosCapabilities-class
### Keywords: classes

### ** Examples

showClass("SosCapabilities_1.0.0")



cleanEx()
nameEx("SosContents")
### * SosContents

flush(stderr()); flush(stdout())

### Name: SosContents-class
### Title: Class and Construction Function of "SosContents"
### Aliases: SosContents-class SosContentsOrNULL-class
###   show,SosContents-method SosContents print,SosContents-method
###   toString,SosContents-method
### Keywords: classes

### ** Examples

showClass("SosContents")
showClass("SosContentsOrNULL")



cleanEx()
nameEx("SosEventTime")
### * SosEventTime

flush(stderr()); flush(stdout())

### Name: SosEventTime
### Title: Classes and Construction Functions for sos:eventTime elements.
### Aliases: SosEventTime-class SosEventTimeLatest-class
###   show,SosEventTime-method print,SosEventTime-method
###   toString,SosEventTime-method SosEventTime SosEventTimeLatest
### Keywords: classes

### ** Examples

showClass("SosEventTime")
showClass("SosEventTimeLatest")

## Not run: 
##D # create SosEventTime for all times after the given time stamp
##D tOps <- TM_After(time = GmlTimeInstant(timePosition = GmlTimePosition(as.POSIXct("2010-01-01 12:00"))))
##D time1 <- SosEventTime(tOps)
##D 
##D # encode it as XML and KVP
##D encodeXML(time1)
##D encodeKVP(time1)
##D 
##D time2 <- SosEventTimeLatest()
##D encodeXML(time2)
## End(Not run)




cleanEx()
nameEx("SosFeatureOfInterest")
### * SosFeatureOfInterest

flush(stderr()); flush(stdout())

### Name: SosFeatureOfInterest-class
### Title: Class and Construction Function for "SosFeatureOfInterest"
### Aliases: SosFeatureOfInterest-class show,SosFeatureOfInterest-method
###   print,SosFeatureOfInterest-method
###   toString,SosFeatureOfInterest-method SosFeatureOfInterest
###   SosFeatureOfInterestOrNULL-class
### Keywords: classes

### ** Examples

showClass("SosFeatureOfInterest")
showClass("SosFeatureOfInterestOrNULL")



cleanEx()
nameEx("SosFilter_Capabilities")
### * SosFilter_Capabilities

flush(stderr()); flush(stdout())

### Name: SosFilter_Capabilities-class
### Title: Classes and Construction Functions for "SosFilter_Capabilities"
###   Elements
### Aliases: SosFilter_Capabilities-class
###   show,SosFilter_Capabilities-method
###   print,SosFilter_Capabilities-method
###   toString,SosFilter_Capabilities-method SosFilter_Capabilities
###   SosFilter_CapabilitiesOrNULL-class
### Keywords: classes

### ** Examples

showClass("SosFilter_Capabilities")
showClass("SosFilter_CapabilitiesOrNULL")



cleanEx()
nameEx("SosObservationOffering")
### * SosObservationOffering

flush(stderr()); flush(stdout())

### Name: SosObservationOffering-class
### Title: Classes and Related Functions for "SosObservationOffering"
### Aliases: SosObservationOffering SosObservationOffering-class
###   show,SosObservationOffering-method
###   sosTime,SosObservationOffering-method
###   sosBoundedBy,SosObservationOffering-method
###   print,SosObservationOffering-method
###   toString,SosObservationOffering-method
###   sosName,SosObservationOffering-method
###   sosResultModels,SosObservationOffering-method
###   sosId,SosObservationOffering-method
### Keywords: classes

### ** Examples

showClass("SosObservationOffering")
# TBD examples for construction functions



cleanEx()
nameEx("Supported")
### * Supported

flush(stderr()); flush(stdout())

### Name: Supported
### Title: Functions to Access Supported Features of the Current sos4R
###   Implementation
### Aliases: SosSupported SosSupportedComparisonOperators
###   SosSupportedConnectionMethods SosSupportedGeometryOperands
###   SosSupportedResponseFormats SosSupportedResponseModes
###   SosSupportedResultModels SosSupportedSpatialOperators
###   SosSupportedTemporalOperators SosSupportedServiceVersions
###   SosSupportedOperations
### Keywords: utilities

### ** Examples

# The supported operations of the specification
SosSupportedOperations()

# HTTP connection methods supported by this sos4R implementation
SosSupportedConnectionMethods()
myConnectionType <- SosSupportedConnectionMethods()[[1]]
myConnectionType

# Formats, modes and models that can be processed by this implementation
SosSupportedResponseFormats()
SosSupportedResultModels()
SosSupportedResponseModes()

# Operators and operands for filtering in a GetObservation request
SosSupportedTemporalOperators()
SosSupportedSpatialOperators()
SosSupportedGeometryOperands()
SosSupportedComparisonOperators()



cleanEx()
nameEx("TM_Operators")
### * TM_Operators

flush(stderr()); flush(stdout())

### Name: TM_Operators
### Title: Classes and Construction Methods for Temporal Operator Classes
### Aliases: print,TM_After-method print,TM_Before-method
###   print,TM_During-method print,TM_Equals-method show,TM_After-method
###   show,TM_Before-method show,TM_During-method show,TM_Equals-method
###   toString,TM_After-method toString,TM_Before-method
###   toString,TM_During-method toString,TM_Equals-method TM_After-class
###   TM_Before-class TM_During-class TM_Equals-class TM_Operators-class
###   TM_Operators TM_After TM_Before TM_During TM_Equals
### Keywords: classes utilities

### ** Examples

showClass("TM_After")
showClass("TM_Before")
showClass("TM_During")
showClass("TM_Equals")

## Not run: 
##D # create times to use for operators
##D t1 <- sosCreateTimeInstant(sos = weathersos, time = Sys.time())
##D p1 <- sosCreateTimePeriod(sos = weathersos, begin = as.POSIXct("2010-03-01 12:15"), end = as.POSIXct("2010-03-02 12:15"))
##D 
##D # create temporal operator
##D afterNow <- TM_After(time = t1)
##D print(afterNow)
##D encodeXML(t1, sos)
##D 
##D during <- TM_During(time = p1)
##D print(during)
## End(Not run)




cleanEx()
nameEx("getObservation-methods")
### * getObservation-methods

flush(stderr()); flush(stdout())

### Name: getObservation-methods
### Title: Request Observations
### Aliases: getObservation getObservation-methods
###   getObservation,SOS,character-method
###   getObservation,SOS,SosObservationOffering-method
###   getObservation,SOS_1.0.0,character-method
###   getObservation,SOS_1.0.0,SosObservationOffering-method
###   getObservationById getObservationById-methods
###   getObservationById,SOS,character-method
###   getObservationById,SOS_1.0.0,character-method
### Keywords: methods

### ** Examples

## Not run: 
##D # by identifier
##D sos <- SOS(...)
##D os.offerings <- sosOfferings(sos)
##D obsId <- getObservationById(sos = sos, observationId = "o_3508493")
##D 
##D # TODO
##D 
## End(Not run)



cleanEx()
nameEx("parse")
### * parse

flush(stderr()); flush(stdout())

### Name: parse
### Title: Parsing Functions for XML Documents and Elements
### Aliases: parse parseCategoryObservation parseComplexObservation
###   parseComponent parseCompositePhenomenon parseCountObservation
###   parseDataArray parseElementType parseEncoding parseFeatureCollection
###   parseField parseFOI parseGeometryObservation parseMeasure
###   parseMeasurement parseObservation parseObservationCollection parseOM
###   parseOwsException parseOwsExceptionReport parseOwsOperation
###   parseOwsRange parseOwsServiceIdentification parseOwsServiceProvider
###   parsePhenomenonProperty parsePoint parsePosition parseResult
###   parseSamplingPoint parseSamplingTime parseSensorML
###   parseSosCapabilities parseSosFilter_Capabilities
###   parseSosObservationOffering parseTemporalObservation parseTextBlock
###   parseTimeGeometricPrimitiveFromParent parseTimeInstant
###   parseTimeInstantProperty parseTimePeriod parseTimePosition
###   parseTruthObservation parseValues parseObservationProperty sosParse
###   sosParse-methods sosParse,SOS_1.0.0,character-method
###   sosParse,SOS_1.0.0,character,logical-method parseNoParsing parseCSV
###   parseFile parseFile-method parseFile,SOS_1.0.0,character-method
###   parseCoordinate parseCoordinate-method parseVector parseVector-method
###   parseLocation parseLocation-method parseSwePosition
###   parseSwePosition-method
### Keywords: methods misc

### ** Examples

# parsing a XML string to an exception report object
er.doc <- xmlParseDoc("<ows:ExceptionReport xmlns:ows=\"http://www.opengis.net/ows/1.1\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\" version=\"1.0.0\" xsi:schemaLocation=\"http://schemas.opengis.net/ows/1.1.0/owsExceptionReport.xsd\"><ows:Exception exceptionCode=\"VersionNegotiationFailed\" locator=\"AcceptVersions\"><ows:ExceptionText>The parameter 'AcceptVersions' does not contain the version of this SOS: '1.0.0'</ows:ExceptionText></ows:Exception></ows:ExceptionReport>")
er.parsed <- parseOwsExceptionReport(er.doc)
print(er.parsed)
str(er.parsed)

## Not run: 
##D # save and re-parse an observation from file
##D obsId <- getObservationById(sos = mySOS, observationId = "o_3508493",
##D 		saveOriginal = TRUE)
##D .files <- list.files(getwd())
##D .startWithO_ <- .files %in% grep("o_", .files, value=TRUE)
##D .observationFiles <- subset(.files, .startWithO_)
##D 
##D obsId <- parseFile(sos = mySOS, file = .observationFiles[[1]])
## End(Not run)




cleanEx()
nameEx("read.sos")
### * read.sos

flush(stderr()); flush(stdout())

### Name: read.sos
### Title: Read Data from a SOS Connection
### Aliases: read.sos
### Keywords: utilities

### ** Examples

# TBD



cleanEx()
nameEx("sos4R-package")
### * sos4R-package

flush(stderr()); flush(stdout())

### Name: sos4R-package
### Title: A client for the OGC Sensor Observation Service
### Aliases: sos4R-package sos4R sosChanges sosCheatSheet sosNews
### Keywords: package connection ts spatial database

### ** Examples


## Not run: 
##D 
##D # Take a SOS from the example list
##D sos.url = SosExampleServices()[[1]]
##D 
##D # Open the connection
##D sos = SOS(url = SOS)
##D 
##D # List offerings, procedures and observedProperties
##D names(sosOfferings(sos))
##D sosProcedures(sos)
##D sosObservedProperties(sos)
##D 
##D # Create time period (last 30 days)
##D tPeriod <- sosCreateEventTimeList(
##D 	time = sosCreateTimePeriod(
##D 		sos = pegelsos,
##D 		begin = Sys.time() - (3600 * 24 * 30),
##D 		end = Sys.time()))
##D 
##D # Request data for all observed properties and procedures of a certain offering
##D observation <- getObservation(sos = sos,
##D 		observedProperty = sosObservedProperties(sos),
##D 		offering = sosOfferings(sos)[[2]],
##D 		procedure = sosProcedures(sos),
##D 		eventTime = tPeriod)
##D 
##D # Inspect result
##D sosResult(observation)
##D str(sosResult(observation))
##D 
##D # Inspect attributes of the data fields
##D if(is.list(sosResult(observation))) {
##D 	attributes(sosResult(observation)[,1])
##D }
##D else {
##D 	attributes(sosResult(pegelObs)[,1])
##D }
##D 
##D # Use custom converting function and connection method. This mechanism works the same for encoders and decoders.
##D myConverters <- SosDataFieldConvertingFunctions(
##D 	"myNumericUnit" = sosConvertDouble)
##D mySos <- SOS(sos.url, method = "GET", dataFieldConverters = myConverters)
##D sosDataFieldConverters(mySos)
##D 
##D # get the cheat sheet
##D sosCheatSheet()
##D 
##D # view the NEWS file
##D sosNews()
##D # DEPRECATED: the changes document
##D #sosChanges()
##D 
## End(Not run)




cleanEx()
nameEx("sosConvert")
### * sosConvert

flush(stderr()); flush(stdout())

### Name: sosConvertString
### Title: SOS Conversion functions for Observation Results
### Aliases: sosConvertString sosConvertDouble sosConvertTime
###   sosConvertLogical
### Keywords: utilities

### ** Examples


## Not run: 
##D sos <- SOS(url = SosExampleServices()[[2]])
##D one <- sosConvertDouble("1", sos)
##D class(one)
##D 
##D # add conversion rules, also possible to override default ones
##D myConverters <- SosDataFieldConvertingFunctions(
##D 	"C" = sosConvertDouble,
##D 	"S/m" = sosConvertDouble)
##D sos <- SOS(url = SosExampleServices()[[2]], dataFieldConverters = myConverters)
##D 
##D # show converters
##D sosDataFieldConverters(sos)
## End(Not run)




cleanEx()
nameEx("sosCreate")
### * sosCreate

flush(stderr()); flush(stdout())

### Name: sosCreate
### Title: Convenience Functions for Request Parameter Creations
### Aliases: sosCreate sosCreateBBOX
###   sosCreateBBOX,numeric,numeric,numeric,numeric-method
###   sosCreateBBoxMatrix
###   sosCreateBBoxMatrix,numeric,numeric,numeric,numeric-method
###   sosCreateEventTimeList sosCreateEventTimeList-methods
###   sosCreateEventTimeList,GmlTimeGeometricPrimitive-method
###   sosCreateEventTime sosCreateEventTime-methods
###   sosCreateEventTime,GmlTimeGeometricPrimitive-method
###   sosCreateFeatureOfInterest sosCreateFeatureOfInterest-methods
###   sosCreateFeatureOfInterest,ANY-method sosCreateTimeInstant
###   sosCreateTimeInstant-methods sosCreateTimeInstant,SOS,POSIXt-method
###   sosCreateTimePeriod sosCreateTimePeriod-methods
###   sosCreateTimePeriod,SOS,POSIXt,POSIXt-method sosCreateTime
###   sosCreateTime-methods sosCreateTime,SOS,character-method
### Keywords: utilities methods

### ** Examples


# create a feature of interest based on identifiers
foiIDs <- list("urn:ogc:object:feature:1", "urn:ogc:object:feature:2")
foiObj <- sosCreateFeatureOfInterest(objectIDs = foiIDs[1:2])
print(foiObj)

# create a bounding box matrix and use it to create a spatial feature of interest
bboxMatrix <- sosCreateBBoxMatrix(lowLat = 50.0, lowLon = 7.0, uppLat = 53.0, uppLon = 10.0)
foiBBox <- sosCreateFeatureOfInterest(bbox = bboxMatrix, srsName = "urn:ogc:def:crs:EPSG:6.8:4326")
print(foiBBox)

# create a foi with a bounding box
bbox <- sosCreateBBOX(lowLat = 50.0, lowLon = 7.0, uppLat = 53.0, uppLon = 10.0, srsName = "urn:ogc:def:crs:EPSG:6.8:4326", srsDimension = as.integer(2), axisLabels = "lat,lon", uomLabels = "deg,deg", propertyName = "bboxName")
foiBBox2 <- sosCreateFeatureOfInterest(spatialOps = bbox)
print(foiBBox2)

## Not run: 
##D last.period <- sosCreateTimePeriod(sos = mySOS,
##D 	begin = (Sys.time() - 3600 * 24 * 7), end = Sys.time())
##D 
##D oneWeek.period <- sosCreateTimePeriod(sos = mySOS,
##D 		begin = as.POSIXct("2010/01/01"), end = as.POSIXct("2010/01/07"))
##D oneWeek.eventTime <- sosCreateEventTimeList(oneWeek.period)
##D 
##D sosCreateTime(sos = mySOS, time = "2007-07-07 07:00::2008-08-08 08:00")
##D sosCreateTime(sos = mySOS, time = "2007-07-07 07:00/2010-10-10 10:00")
##D 
##D sosCreateTime(sos = mySOS, time = "::2007-08-05")
##D sosCreateTime(sos = mySOS, time = "2007-08-05/")
## End(Not run)




### * <FOOTER>
###
base::cat("Time elapsed: ", proc.time() - base::get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
