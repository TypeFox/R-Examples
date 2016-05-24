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
# Created: 2010-06-18                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#
# Constants for version 1.0.0 of SOS
#

################################################################################
# SOS
sosService <- "SOS"
sosNamespacePrefix <- "sos"

# Core Operations Profile:
sosGetCapabilitiesName <- "GetCapabilities"
sosDescribeSensorName <- "DescribeSensor"
sosGetObservationName <- "GetObservation"
# Transaction Operations Profile
sosRegisterSensorName <- "RegisterSensor"
sosInsertObservationName <- "InsertObservation"
# Enhanced Operations Profile:
sosGetObservationByIdName <- "GetObservationById"
sosGetResultName <- "GetResult"
sosGetFeatureOfInterestName <- "GetFeatureOfInterest"
sosGetFeatureOfInterestTimeName <- "GetFeatureOfInterestTime"
sosDescribeFeatureTypeName <- "DescribeFeatureType"
sosDescribeObservationTypeName <- "DescribeObservationType"
sosDescribeResultModelName <- "DescribeResultModel"

SosSupportedOperations <- function() {
	.supported <- c(sosGetCapabilitiesName, sosDescribeSensorName, 
			sosGetObservationName ,sosGetObservationByIdName)
	return(.supported)
}

################################################################################
# not exported SOS
.sosConnectionMethodGet <- "GET"
.sosConnectionMethodPost <- "POST"
.sosConnectionMethodSOAP <- "SOAP"

SosSupportedConnectionMethods <- function() {
	.supported <- c(.sosConnectionMethodGet, .sosConnectionMethodPost)
	names(.supported) <- c(.sosConnectionMethodGet, .sosConnectionMethodPost)
	return(.supported)
}

mimeTypeCSV <- "text/csv"
mimeTypeXML <- "text/xml"
mimeTypeOM <- "text/xml;subtype=&quot;om/1.0.0&quot;"
mimeTypeSML <- "text/xml;subtype=&quot;sensorML/1.0.1&quot;"
mimeTypeKML <- "application/vnd.google-earth.kml+xml"
mimeSubtypeOM <- "\"om/1.0.0\""

.sosSupportedResponseFormats <- c(
		mimeTypeOM,
		mimeTypeSML,
		mimeTypeCSV,
		mimeTypeKML)
SosSupportedResponseFormats <- function() {
	return(.sosSupportedResponseFormats)
}

.sosSupportedResultModels <- c("om:Measurement", "om:Observation")
SosSupportedResultModels <- function() {
	return(.sosSupportedResultModels)
}

.sosSupportedResponseModes <- c("inline")
SosSupportedResponseModes <- function() {
	return(.sosSupportedResponseModes)
}

.sosSupportedServiceVersions <- c("1.0.0")
SosSupportedServiceVersions <- function() {
	return(.sosSupportedServiceVersions)
}

.sosNamespaceDefinitionsForAll <- c(sos = "http://www.opengis.net/sos/1.0",
		xsi = "http://www.w3.org/2001/XMLSchema-instance")
.sosNamespaceDefinitionsGetObs <- c(ows = "http://www.opengis.net/ows/1.1",
		om = "http://www.opengis.net/om/1.0",
		ogc = "http://www.opengis.net/ogc",
		gml = "http://www.opengis.net/gml")
.sosNamespaceDefinitionsGetCap <- c(ows = "http://www.opengis.net/ows/1.1",
		ogc = "http://www.opengis.net/ogc")
.sosNamespaceDefinitionsSML <- c(sml = "http://www.opengis.net/sensorML/1.0.1",
		gml = "http://www.opengis.net/gml",
		swe = "http://www.opengis.net/swe/1.0.1",
		xlink = "http://www.w3.org/1999/xlink",
		xsi = "http://www.w3.org/2001/XMLSchema-instance")
		

.xsiSchemaLocationAttribute <- c("xsi:schemaLocation" = "http://www.opengis.net/sos/1.0 http://schemas.opengis.net/sos/1.0.0/sosAll.xsd")

################################################################################
# SOS
sosIntendedApplicationName <- "intendedApplication"
sosTimeName <- "time"
sosProcedureName <- "procedure"
sosObservedPropertyName <- "observedProperty"
sosFeatureOfInterestName <- "featureOfInterest"
sosResultModelName <- "resultModel"
sosResponseFormatName <- "responseFormat"
sosResponseModeName <- "responseMode"
sosObservationOfferingName <- "ObservationOffering"
sosObservationOfferingListName <- "ObservationOfferingList"
sosContentsName <- "Contents"
sosFilterCapabilitiesName <- "Filter_Capabilities"
sosCapabilitiesName <- "Capabilities"
sosEventTimeName <- "eventTime"
sosEventTimeLatestValue <- "latest"
sosObjectIDName <- "ObjectID"
sosResultName <- "result"

################################################################################
# O&M
omMeasurementName <- "Measurement"
omMemberName <- "member"
omObservationName <- "Observation"
omObservationCollectionName <- "ObservationCollection"
omFeatureOfInterestName <- "featureOfInterest"
omProcedureName <- "procedure"
omObservedPropertyName <- "observedProperty"
omSamplingTimeName <- "samplingTime"
omResultTimeName <- "resultTime"
omResultName <- "result"
omCategoryObservationName <- "CategoryObservation"
omCountObservationName <- "CountObservation"
omTruthObservationName <- "TruthObservation"
omGeometryObservationName <- "GeometryObservation"
omTemporalObservationName <- "TemporalObservation"
omComplexObservationName <- "ComplexObservation"

################################################################################
# SA
saSamplingPointName <- "SamplingPoint"
saSamplingSurface <- "SamplingSurface"
saPositionName <- "position"
saSampledFeatureName <- "sampledFeature"
saSamplingTimeName <- "samplingTime"

################################################################################
# GML
gmlPosName <- "pos"
gmlPointName <- "Point"
gmlTimeInstantName <- "TimeInstant"
gmlTimePeriodName <- "TimePeriod"
gmlTimePositionName <- "timePosition"
gmlRelatedTimeName <- "relatedTime"
gmlNameName <- "name"
gmlDescriptionName <- "description"
gmlBeginName <- "begin"
gmlEndName <- "end"
gmlBeginPositionName <- "beginPosition"
gmlEndPositionName <-"endPosition"
gmlFeatureCollectionName <- "FeatureCollection"
gmlBoundedByName <- "boundedBy"
gmlEnvelopeName <- "Envelope"
gmlLowerCornerName <- "lowerCorner"
gmlUpperCornerName <- "upperCorner"
gmlNamespacePrefix <- "gml"
gmlTimeLengthName <- "timeLength"
gmlDurationName <- "duration"
gmlTimeIntervalName <- "timeInterval"
gmlFeatureMemberName <- "featureMember"

################################################################################
# SWE
sweCompositePhenomenonName <- "CompositePhenomenon"
sweBaseName <- "base"
sweComponentName <- "component"
sweDataArrayName <- "DataArray"
sweElementTypeName <- "elementType"
sweSimpleDataRecordName <- "SimpleDataRecord"
sweDataRecordName <- "DataRecord"
sweFieldName <- "field"
sweTimeName <- "Time"
sweQuantityName <- "Quantity"
sweCategoryName <- "Category"
sweBooleanName <- "Boolean"
sweCountName <- "Count"
sweEncodingName <- "encoding"
sweTextBlockName <- "TextBlock"
sweValuesName <- "values"
sweValueName <- "value"
sweCodeSpaceName <- "codeSpace"
sweTextName <- "Text"
sweUomName <- "uom"
sweVectorName <- "Vector"
sweLocationName <- "location"
sweCoordinateName <- "coordinate"
swePositionName <- "Position"

################################################################################
# OGC
ogcTempOpTMAfterName <- "TM_After"
ogcTempOpTMBeforeName <- "TM_Before"
ogcTempOpTMBeginsName <- "TM_Begins"
ogcTempOpTMBegunByName <- "TM_BegunBy"
ogcTempOpTMContainsName <- "TM_Contains"
ogcTempOpTMDuringName <- "TM_During"
ogcTempOpTMEndedByName <- "TM_EndedBy"
ogcTempOpTMEndsName <- "TM_Ends"
ogcTempOpTMEqualsName <- "TM_Equals"
ogcTempOpTMMeetsName <- "TM_Meets"
ogcTempOpTMMetByName <- "TM_MetBy"
ogcTempOpTMOverlapsName <- "TM_Overalps"
ogcTempOpTMOverlappedBy <- "TM_OverlappedBy"
.ogcSupportedTemporalOps <- list(
		ogcTempOpTMAfterName,
		ogcTempOpTMBeforeName,
		ogcTempOpTMDuringName,
		ogcTempOpTMEqualsName
)
names(.ogcSupportedTemporalOps) <- .ogcSupportedTemporalOps
SosSupportedTemporalOperators <- function() {
	return(.ogcSupportedTemporalOps)
}

ogcSpatialOpBBOXName <- "BBOX"
ogcSpatialOpContainsName <- "Contains"
ogcSpatialOpIntersectsName <- "Intersects"
ogcSpatialOpOverlapsName <- "Overlaps"
ogcSpatialOpBeyondName <- "Beyond"
ogcSpatialOpCrossesName <- "Crosses"
ogcSpatialOpDWithinName <- "DWithin"
ogcSpatialOpDisjointName <- "Disjoint"
ogcSpatialOpEqualsName <- "Equals"
ogcSpatialOpTouchesName <- "Touches"
ogcSpatialOpWithinName <- "Within"
.ogcSupportedSpatialOps <- list(
		ogcSpatialOpBBOXName,
		ogcSpatialOpContainsName,
		ogcSpatialOpIntersectsName,
		ogcSpatialOpOverlapsName
)
names(.ogcSupportedSpatialOps) <- .ogcSupportedSpatialOps
SosSupportedSpatialOperators <- function() {
	return(.ogcSupportedSpatialOps)
}

ogcGeometryOperandEnvelopeName <- "gml:Envelope"
ogcGeometryOperandPolygonName <- "gml:Polygon"
ogcGeometryOperandPointName <- "gml:Point"
ogcGeometryOperandLineStringName <- "gml:LineString"

.ogcSupportedGeometryOperands <- list(
		ogcGeometryOperandEnvelopeName,
		ogcGeometryOperandPolygonName,
		ogcGeometryOperandPointName,
		ogcGeometryOperandLineStringName
)
names(.ogcSupportedGeometryOperands) <- .ogcSupportedGeometryOperands
SosSupportedGeometryOperands <- function() {
	return(.ogcSupportedGeometryOperands)
}

ogcComparisonOpBetweenName <- "PropertyIsBetween"
ogcComparisonOpEqualToName <- "PropertyIsEqualTo"
ogcComparisonOpGreaterThanName <- "PropertyIsGreaterThan"
ogcComparisonOpGreaterThanOrEqualToName <- "PropertyIsGreaterThanOrEqualTo"
ogcComparisonOpLessThenName <- "PropertyIsLessThan"
ogcComparisonOpLessThanOrEqualToName <- "PropertyIsLessThanOrEqualTo"
ogcComparisonOpIsLikeName <- "PropertyIsLike"
ogcComparisonOpIsNotEqualTo <- "PropertyIsNotEqualTo"
ogcComparisonOpIsNull <- "PropertyIsNull"
.ogcSupportedComparisonOperators <- list()
names(.ogcSupportedComparisonOperators) <- .ogcSupportedComparisonOperators
SosSupportedComparisonOperators <- function() {
	return(.ogcSupportedComparisonOperators)
}

ogcNamespacePrefix <- "ogc"
ogcPropertyNameName <- "PropertyName"
ogcBBOXName <- "BBOX"
ogcContainsName <- "Contains"
ogcIntersectsName <- "Intersects"
ogcOverlapsName <- "Overlaps"
ogcSpatialCapabilitiesName <- "Spatial_Capabilities"
ogcTemporalCapabilitiesName <- "Temporal_Capabilities"
ogcScalarCapabilitiesName <- "Scalar_Capabilities"
ogcIdCapabilities <- "Id_Capabilities"
ogcGeometryOperandsName <- "GeometryOperands"
ogcGeometryOperandName <- "GeometryOperand"
ogcSpatialOperatorsName <- "SpatialOperators"
ogcSpatialOperatorName <- "SpatialOperator"
ogcTemporalOperandsName <- "TemporalOperands"
ogcTemporalOperandName <- "TemporalOperand" 
ogcTemporalOperatorsName <- "TemporalOperators"
ogcTemporalOperatorName <- "TemporalOperator"
ogcLogicalOperatorsName <- "LogicalOperators"
ogcComparisonOperatorsName <- "ComparisonOperators"
ogcArithmeticOperatorsName <- "ArithmeticOperators"
ogcEIDName <- "EID"
ogcFIDName <- "FID"
ogcLiteralName <- "Literal"

smlSensorMLName <- "SensorML"

################################################################################
# OWS
owsServiceIdentificationName <- "ServiceIdentification"
owsTitleName <- "Title"
owsAbstractName <- "Abstract"
owsKeywordsName <- "Keywords"
owsKeywordName <- "Keyword"
owsServiceTypeName <- "ServiceType"
owsServiceTypeVersionName <- "ServiceTypeVersion"
owsFeesName <- "Fees"
owsAccessConstraintsName <- "AccessConstraints"
owsServiceProviderName <- "ServiceProvider"
owsOperationsMetadataName <- "OperationsMetadata"
owsOperationName <- "Operation"
owsDCPName <- "DCP"
owsHTTPName <- "HTTP"
owsGetName <- "Get"
owsPostName <- "Post"
owsParameterName <- "Parameter"
owsAllowedValuesName <- "AllowedValues"
owsValueName <- "Value"
owsAnyValueName <- "AnyValue"
owsRangeName <- "Range"
owsMinimumValueName <- "MinimumValue"
owsMaximumValueName <- "MaximumValue"
owsSpacingName <- "Spacing"
owsConstraintName <- "Constraint"
owsMetadataName <- "Metadata"
owsExceptionReportName <- "ExceptionReport"
owsExceptionName <- "Exception"
owsExceptionTextName <- "ExceptionText"
owsProfileName <- "Profile"
owsProviderNameName <- "ProviderName"
owsProviderSiteName <- "ProviderSite"
owsServiceContactName <- "ServiceContact"

kmlName <- "kml"

################################################################################
owsNamespacePrefix <- "ows"
.owsNamespace <- c(ows = "http://www.opengis.net/ows/1.1")
.owsCodes = c(
		"OperationNotSupported",
		"MissingParameterValue",
		"InvalidParameterValue",
		"VersionNegotiationFailed",
		"InvalidUpdateSequence",
		"OptionNotSupported",
		"NoApplicableCode")
.owsCodeMeanings = c(
		"Request is for an operation that is not supported by this server",
		"Operation request does not include a parameter value, and this server did not declare a default parameter value for that parameter",
		"Operation request contains an invalid parameter value",
		"List of versions in 'AcceptVersions' parameter value in GetCapabilities operation request did not include any version supported by this server",
		"Value of (optional) updateSequence parameter in GetCapabilities operation request is greater than current value of service metadata updateSequence number",
		"Request is for an option that is not supported by this server",
		"No other exceptionCode specified by this service and server applies to this exception")
.owsCodeLocators = c(
		"Name of operation not supported",
		"Name of missing parameter",
		"Name of parameter with invalid value",
		"None, omit 'locator' parameter",
		"None, omit 'locator' parameter",
		"Identifier of option not supported",
		"None, omit 'locator' parameter")
.httpCode = c("501", "400", "400", "400", "400", "501", "3xx, 4xx, 5xx")
.httpMessage = c("Not Implemented", "Bad request", "Bad request", "Bad request",
		"Bad request", "Not implemented", "Internal Server Error")

.owsStandardExceptions <- data.frame(
		exceptionCode = .owsCodes,
		meaningOfCode = .owsCodeMeanings, 
		locator = .owsCodeLocators,
		httpStatusCode = .httpCode,
		httpMessage = .httpMessage,
		check.rows = TRUE, check.names = TRUE)
OwsExceptionsData <- function() {
	return(.owsStandardExceptions)
}

################################################################################
# others
xmlInternalDocumentName <- "XMLInternalDocument"
xmlTextNodeName <- "text"

.sosCheatSheetDocumentName <- "sos4r_cheat-sheet.pdf"
sosAttributeFileName <- "savedAsFile"
