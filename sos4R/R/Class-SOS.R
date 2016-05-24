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
#
#
setClass("SOS",
		representation(version = "character",
				capabilities = "OwsCapabilities", parsers = "list",
				encoders = "list", dataFieldConverters = "list",
				timeFormat = "character", verboseOutput = "logical",
				switchCoordinates = "logical"),
		contains = c("VIRTUAL"))

#
#
#
setClass("SOS_1.0.0",
		representation(url = "character", method = "character",
				curlHandle = "CURLHandle", curlOptions = "ANY"),
		prototype = list(
				url = as.character(NA),
				method = as.character(NA),
				version = as.character(NA)),
		contains = c("SOS"),
		validity = function(object) {
			#print("Entering validation: SOS")
			
			# method has to be one of GET, POST, SOAP
			.allowedMethods = c(.sosConnectionMethodGet,
					.sosConnectionMethodPost, .sosConnectionMethodSOAP)
			if(!any(sapply(.allowedMethods, "==", object@method), na.rm = TRUE))
				return(paste("method has to be one of", paste(.allowedMethods,
										sep = ", ", collapse = " ")))
				
			# url has to match an URL pattern
			.urlPattern = "(?:https?://(?:(?:(?:(?:(?:[a-zA-Z\\d](?:(?:[a-zA-Z\\d]|-)*[a-zA-Z\\d])?)\\.)*(?:[a-zA-Z](?:(?:[a-zA-Z\\d]|-)*[a-zA-Z\\d])?))|(?:(?:\\d+)(?:\\.(?:\\d+)){3}))(?::(?:\\d+))?)(?:/(?:(?:(?:(?:[a-zA-Z\\d$\\-_.+!*'(),]|(?:%[a-fA-F\\d]{2}))|[;:@&=])*)(?:/(?:(?:(?:[a-zA-Z\\d$\\-_.+!*'(),]|(?:%[a-fA-F\\d]{2}))|[;:@&=])*))*)(?:\\?(?:(?:(?:[a-zA-Z\\d$\\-_.+!*'(),]|(?:%[a-fA-F\\d]{2}))|[;:@&=])*))?)?)"
			.result = regexpr(.urlPattern, object@url)
			if (.result == -1)
				return("url not matching URL-pattern (http://www.example.com)")
			
			# test for complete match removed, does not work yet
			#.urlLength = nchar(object@url)
			#if (.urlLength == attr(.result, "match.length"))
			#	return("url not completely matching URL-pattern")
				
			return(TRUE)
		}
)

#
#
#
setClass("SosFilter_Capabilities",
		representation(spatial = "list", temporal = "list", scalar = "list",
				id = "list"),
		prototype = list(spatial = list(NA_character_),
				temporal = list(NA_character_), scalar = list(NA_character_),
				id = list(NA_character_)),
		validity = function(object) {
			#print("Entering validation: SosFilter_Capabilities")
			# TODO implement validity function
			return(TRUE)
		}
)
setClassUnion(name = "SosFilter_CapabilitiesOrNULL",
		members = c("SosFilter_Capabilities", "NULL"))

#
# See OWS Common 1.1.0, OGC 06-121r3
#
setClass("SosCapabilities_1.0.0",
		representation(filterCapabilities = "SosFilter_CapabilitiesOrNULL"),
		contains = "OwsCapabilities_1.1.0",
		validity = function(object) {
			#print("Entering validation: SosCapabilities_1.0.0")
			# TODO implement validity function
			return(TRUE)
		}
)

#
# See OGC 06-009r6, clause 8.2.3.2
#
# Intermediate element swe:TimeGeometricPrimitiveProperty for time is left out.
#
# If possible, character vectors are used for slots if elements can be
# represented by character strings, e.g. procedure is gml:ReferenceType,
# resultModel is QName, intendedApplication is xs:token, and so forth.
#
setClass("SosObservationOffering",
		representation(id = "character", name = "character",
				time = "GmlTimeGeometricPrimitive", procedure = "character",
				observedProperty = "list", featureOfInterest = "list",
				responseFormat = "character", intendedApplication = "character",
				resultModel = "character", responseMode = "character",
				boundedBy = "list"),
		prototype = list(id = as.character(NA), time = NULL,
				procedure = as.character(NA), observedProperty = list(NA),
				featureOfInterest = list(NA),
				responseFormat =  as.character(NA)),
		validity = function(object) {
			#print("Entering validation: ObservationOffering")
			# TODO implement validity function
			
			# time is required
			# procedure, observedProperty, featureOfInterest, responseFormat are all "one or more"
			
			return(TRUE)
		}
)

#
# See OGC 06-009r6, clause 8.2.3.2
#
setClass("SosContents",
		representation(observationOfferings = "list"),
		prototype = list(observationOfferings = list(NA), xml = xmlNode(NA)),
		contains = c("OwsContents"),
		validity = function(object) {
			#print("Entering validation: SosContents")
			# TODO implement validity function
			return(TRUE)
		}
)
setClassUnion(name = "SosContentsOrNULL", members = c("SosContents", "NULL"))

#
#
#
setClass("SosEventTime",
		representation(temporalOps = "OgcBinaryTemporalOp"),
		prototype = list(temporalOps = NULL),
		validity = function(object) {
			#print("Entering validation: SosEventTime")
			# TODO implement validity function
			return(TRUE)
		}
)

#
# class needed for the non-standard request for latest value
#
setClass("SosEventTimeLatest",
		representation(temporalOps = "character"),
		prototype = list(temporalOps = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: SosEventTimeLatest")
			# TODO implement validity function
			return(TRUE)
		}
)

#
# 
#
setClass("SosFeatureOfInterest",
		representation(
				objectIDs = "list",
				spatialOps = "OgcSpatialOpsOrNULL"),
		prototype = list(
				objectIDs = list(NA),
				spatialOps = NULL),
		validity = function(object) {
			#print("Entering validation: GetObservationById")
			# TODO implement validity function
			# one of objectIDs or spatialOps has to be set
			return(TRUE)
		}
)
setClassUnion(name = "SosFeatureOfInterestOrNULL",
		members = c("SosFeatureOfInterest", "NULL"))


################################################################################
# See OGC 06-009r6
#
setClass("SosDescribeSensor",
		representation(procedure = "character", outputFormat = "character"),
		prototype = list(service = as.character(NA), version = as.character(NA),
				procedure = as.character(NA), outputFormat = as.character(NA)),
		contains = "OwsServiceOperation",
		validity = function(object) {
			#print("Entering validation: sosDescribeSensor")
			# TODO implement validity function
			# check format of version, sensorid and outputformat?!
			if(length(object@procedure) > 1)
				return("can only request one procedure at a time!")
			
			return(TRUE)
		}
)


################################################################################
# See SOS specification, OGC 06-009r6, section 8.4
#
setClass("SosGetObservation",
		representation(
				offering = "character",
				observedProperty = "list",
				responseFormat = "character",
				srsName = "character",
				eventTime = "list", 
				procedure = "character", 
				featureOfInterest = "SosFeatureOfInterestOrNULL", 
				result = "ANY", # OgcComparisonOpsOrXMLOrNULL
				resultModel = "character",
				responseMode = "character",
				BBOX = "character"),
		prototype = list(
				service = as.character(NA),
				version = as.character(NA),
				offering = as.character(NA),
				observedProperty = list(NA),
				responseFormat = as.character(NA)),
		contains = "OwsServiceOperation",
		validity = function(object) {
			#print("Entering validation: SosGetObservation")
			# TODO implement validity function
			
			# service, version, offering, observedProperty, and responseFormat are mandatory
			if(is.na(object@service))
				return("service parameter must be given")
			if(is.na(object@version))
				return("version must be given")
			if(is.na(object@offering))
				return("offering parameter must be given")
			# responseFormat is optional for GET
			#if(is.na(object@responseFormat))
			#	return("responseFormat parameter must be given")
			if(length(object@observedProperty) < 1)
				return("at least one observedProperty is mandatory")
			
			# if version is there, it hast to be in a certain format, see ows common
			# srsName, offering, procedure, observedProperty are anyURIs
			# eventTime is a list of ogc:temporalOps
			# featureOfInterest is null or a SosFeatureOfInterest element
	
			# result is null or an ogc:comparisonOps element
			cls <- class(slot(object, "result"))
#			print(paste("class of result slot: ", cls))
			if ( !any(cls %in% c("OgcComparisonOps", "XMLNode", "NULL",
							"XMLAbstractNode", "XMLInternalNode"))) {
				return("'response' argument does not have allowed class!")
			}
			
			# responseFormat must be MIME content type
			# resultModel must be a QName
			# responseMode must be one of inline, out-of-band, attached, or resultTemplate
			return(TRUE)
		}
)


################################################################################
# See SOS specification, OGC 06-009r6, section 10.1
#
setClass("SosGetObservationById",
		representation(
				observationId = "character",
				responseFormat = "character",
				srsName = "character",
				resultModel = "character",
				responseMode = "character"),
		prototype = list(
				observationId = as.character(NA),
				responseFormat = as.character(NA)),
		contains = "OwsServiceOperation",
		validity = function(object) {
			#print("Entering validation: SosGetObservationById")
			# TODO implement validity function
			
			# service, version, observationId, and responseFormat are mandatory
			if(is.na(object@service))
				return("service parameter must be given")
			if(is.na(object@version))
				return("version must be given")
			if(is.na(object@observationId))
				return("observationId parameter must be given")
			if(is.na(object@responseFormat))
				return("responseFormat parameter must be given")
			
			# if version is there, it hast to be in a certain format, see ows common
			# responseFormat must be MIME content type
			# resultModel must be a QName
			# responseMode must be one of inline, out-of-band, attached, or resultTemplate
			return(TRUE)
		}
)
