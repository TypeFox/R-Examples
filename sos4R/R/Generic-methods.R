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
# visit the Free Software Foundation web page, http://www.fsf.org              #
#                                                                              #
# Author: Daniel Nuest (daniel.nuest@uni-muenster.de)                          #
# Created: 2010-09-20                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################


#
#
#
if (!isGeneric("sosRequest"))
	setGeneric(name = "sosRequest",
			def = function(sos, request, verbose = sos@verboseOutput,
					inspect = FALSE) {
				standardGeneric("sosRequest")
			})

#
#
#
if (!isGeneric("getCapabilities"))
	setGeneric(name = "getCapabilities",
			signature = signature("sos", "verbose", "inspect"),
			def = function(sos, verbose = sos@verboseOutput, inspect = FALSE,
					sections = sosDefaultGetCapSections,
					acceptFormats = sosDefaultGetCapAcceptFormats,
					updateSequence = c(as.character(NA)),
					owsVersion = sosDefaultGetCapOwsVersion,
					acceptLanguages = c(NA)) {
				standardGeneric("getCapabilities")	
			})

#
#
#
if (!isGeneric("describeSensor"))
	setGeneric(name = "describeSensor",
			signature = signature("sos", "procedure", "outputFormat", "verbose",
					"inspect", "saveOriginal"),
			def = function(sos, procedure,
					outputFormat = sosDefaultDescribeSensorOutputFormat,
					verbose = sos@verboseOutput, inspect = FALSE,
					saveOriginal = NULL) {
				standardGeneric("describeSensor")	
			})

#
#
#
if (!isGeneric("getObservationById"))
	setGeneric(name = "getObservationById",
			signature = signature("sos", "observationId", "responseFormat",
					"srsName", "resultModel", "responseMode", "verbose", 
					"inspect", "saveOriginal"),
			def = function(sos, observationId,
					responseFormat = sosDefaultGetObsResponseFormat,
					srsName = as.character(NA), resultModel = as.character(NA),
					responseMode = as.character(NA),
					verbose = sos@verboseOutput, inspect = FALSE,
					saveOriginal = NULL) {
				standardGeneric("getObservationById")
			})

#
#
#
if (!isGeneric("getObservation"))
	setGeneric(name = "getObservation",
			signature = signature("sos", "offering", "observedProperty",
					"responseFormat", "srsName", "eventTime", "procedure",
					"featureOfInterest", "result", "resultModel",
					"responseMode", "BBOX", "latest", "verbose", "inspect",
					"saveOriginal"),
			def = function(sos, offering,
					observedProperty = sosObservedProperties(obj = offering),
					responseFormat = sosDefaultGetObsResponseFormat,
					# optional:
					srsName = as.character(NA),
					eventTime = list(NA), # sosCreateEventTimeList(time = sosTime(obj = offering))
					procedure = as.character(NA), # sosProcedures(obj = offering),
					featureOfInterest = NULL,
					result = NULL,
					resultModel = as.character(NA),
					responseMode = as.character(NA),
					BBOX = as.character(NA),
					latest = FALSE,
					verbose = sos@verboseOutput,
					inspect = FALSE,
					saveOriginal = NULL) {
				standardGeneric("getObservation")
			})

#
#
#
if (!isGeneric("checkRequest"))
	setGeneric(name = "checkRequest",
			def = function(service, operation, verbose) {
				standardGeneric("checkRequest")
			})

#
#
#
if (!isGeneric("encodeRequestKVP"))
	setGeneric(name = "encodeRequestKVP",
			def = function(obj, sos, verbose = FALSE) {
				standardGeneric("encodeRequestKVP")
			})

#
#
#
if (!isGeneric("encodeRequestXML"))
	setGeneric(name = "encodeRequestXML",
			def = function(obj, sos, verbose = FALSE) {
				standardGeneric("encodeRequestXML")
			})

#
#
#
if (!isGeneric("encodeRequestSOAP"))
	setGeneric(name = "encodeRequestSOAP",
			def = function(obj, sos, verbose = FALSE) {
				standardGeneric("encodeRequestSOAP")
			})

#
#
#
if (!isGeneric("sosExceptionCodeMeaning"))
	setGeneric(name = "sosExceptionCodeMeaning", def = function(exceptionCode) {
				standardGeneric("sosExceptionCodeMeaning")
			})

#
#
#
if (!isGeneric("encodeXML"))
	setGeneric(name = "encodeXML",
			def = function(obj, sos, verbose = FALSE, ...) {
				standardGeneric("encodeXML")
			})

#
#
#
if (!isGeneric("encodeKVP"))
	setGeneric(name = "encodeKVP",
			def = function(obj, sos, verbose = FALSE, ...) {
				standardGeneric("encodeKVP")
			})

#
#
#
if (!isGeneric("sosGetCRS"))
	setGeneric(name = "sosGetCRS",
			def = function(obj, verbose = FALSE) {
				standardGeneric("sosGetCRS")
			})

#
#
#
if (!isGeneric("parseFile"))
	setGeneric(name = "parseFile",
			def = function(sos, file, verbose = FALSE, ...) {
				standardGeneric("parseFile")
			})

#
#
#
if (!isGeneric("sosGetDCP"))
	setGeneric(name = "sosGetDCP",
			def = function(sos, operation, type = NA) {
				standardGeneric("sosGetDCP")
			})

#
#
#
if (!isGeneric("sosCreateEventTime"))
	setGeneric(name = "sosCreateEventTime",
			def = function(time, operator = sosDefaultTemporalOperator) {
				standardGeneric("sosCreateEventTime")
			})

#
#
#
if (!isGeneric("sosCreateTimeInstant"))
	setGeneric(name = "sosCreateTimeInstant", def = function(sos, time,
					frame = as.character(NA),
					calendarEraName = as.character(NA),
					indeterminatePosition = as.character(NA)) {
				standardGeneric("sosCreateTimeInstant")
			}
	)

#
#
#
if (!isGeneric("sosCreateTimePeriod"))
	setGeneric(name = "sosCreateTimePeriod",
			def = function(sos, begin, end, frame = as.character(NA),
					calendarEraName = as.character(NA),
					indeterminatePosition = as.character(NA),
					duration = as.character(NA),
					timeInterval = NULL) {
				standardGeneric("sosCreateTimePeriod")
			}
	)

#
#
#
if (!isGeneric("sosCreateEventTimeList"))
	setGeneric(name = "sosCreateEventTimeList",
			def = function(time, operator = sosDefaultTemporalOperator) {
				standardGeneric("sosCreateEventTimeList")
			})

#
#
#
if (!isGeneric("sosCreateTime"))
	setGeneric(name = "sosCreateTime",
			def = function(sos, time, operator = sosDefaultTemporalOperator) {
				standardGeneric("sosCreateTime")
			})

#
#
#
if (!isGeneric("sosCreateFeatureOfInterest"))
	setGeneric(name = "sosCreateFeatureOfInterest",
			def = function(objectIDs = list(NA), spatialOps = NULL, bbox = NULL,
					srsName = NA_character_) {
				standardGeneric("sosCreateFeatureOfInterest")
			})

#
#
#
if (!isGeneric("sosCreateBBOX"))
	setGeneric(name = "sosCreateBBOX",
			def = function(lowLat, lowLon, uppLat, uppLon, srsName,
					srsDimension = NA_integer_, axisLabels = NA_character_,
					uomLabels = NA_character_,
					propertyName = sosDefaultSpatialOpPropertyName) {
				standardGeneric("sosCreateBBOX")
			})

#
#
#
if (!isGeneric("sosCreateBBoxMatrix"))
	setGeneric(name = "sosCreateBBoxMatrix",
			def = function(lowLat, lowLon, uppLat, uppLon) {
				standardGeneric("sosCreateBBoxMatrix")
			})

#
#
#
if (!isGeneric("sosCapabilitiesDocumentOriginal"))
	setGeneric(name = "sosCapabilitiesDocumentOriginal", def = function(sos) {
				standardGeneric("sosCapabilitiesDocumentOriginal")
			})

#
#
#
if (!isGeneric("sosCapabilitiesUrl"))
	setGeneric(name = "sosCapabilitiesUrl",
			def = function(sos) {
				standardGeneric("sosCapabilitiesUrl")
			})
