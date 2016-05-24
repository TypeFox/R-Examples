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

################################################################################
# construction functions
SOS <- function(url, method = SosDefaultConnectionMethod(),
		version = "1.0.0",
		parsers = SosParsingFunctions(),
		encoders = SosEncodingFunctions(),
		dataFieldConverters = SosDataFieldConvertingFunctions(),
		curlOptions = list(),
		curlHandle = getCurlHandle(),
		timeFormat = sosDefaultTimeFormat, verboseOutput = FALSE, 
		switchCoordinates = FALSE, ...) {
	if(version == "1.0.0") {
		if(method == .sosConnectionMethodPost)
			.curlOpts <- curlOptions(url = url)
		else .curlOpts <- curlOptions
		
		.sos <- new("SOS_1.0.0",
				url = url,
				method = method,
				version = version,
				# dummy capabilities to be replaced below
				capabilities = new("OwsCapabilities", version = "NA",
						updateSequence = as.character(NA),
						owsVersion = sosDefaultGetCapOwsVersion),
				parsers = parsers,
				encoders = encoders,
				dataFieldConverters = dataFieldConverters,
				curlOptions = .curlOpts,
				curlHandle = curlHandle,
				timeFormat = timeFormat,
				verboseOutput = verboseOutput,
				switchCoordinates = switchCoordinates)
		
		.caps <- getCapabilities(sos = .sos, verbose = verboseOutput, ...)
		if(!is(.caps, "OwsCapabilities")) {
			stop("ERROR: Did not receive a Capabilities response!")
		}
		
		.sos@capabilities <- .caps
		
		if(verboseOutput) cat("[SOS] Created new SOS:\n", toString(.sos), "\n")
				
		cat("[sos4R] Created SOS for URL", url, "\n")
		return(.sos)
	}
	else stop("Service version not supported!")
	
}


SosFilter_Capabilities <- function(spatial = list(NA_character_),
		temporal = list(NA_character_), scalar = list(NA_character_),
		id = list(NA_character_)) {
	new("SosFilter_Capabilities", spatial = spatial, temporal = temporal,
			scalar = scalar, id = id)
}

SosCapabilities <- function(version,  updateSequence = NA, owsVersion = "1.1.0",
		identification = NULL, provider = NULL, operations = NULL,
		filterCapabilities = NULL, contents = NULL) {
	if(owsVersion == "1.1.0") {
		new("SosCapabilities_1.0.0",
				version = version, updateSequence = updateSequence,
				owsVersion = owsVersion,
				identification = identification,
				provider = provider, operations = operations,
				filterCapabilities = filterCapabilities, contents = contents)
	}
	else if(owsVersion == "2.0.0") {
		stop("Version 2.0.0 not supported!")
	}
	else {
		new("OwsCapabilities",
				version = version, updateSequence = updateSequence,
				owsVersion = owsVersion)
	}	
}

SosObservationOffering <- function(id, name = as.character(NA),
		time, procedure, observedProperty,
		featureOfInterest, responseFormat,
		intendedApplication = as.character(NA), resultModel = as.character(NA),
		responseMode = as.character(NA), boundedBy = list()) {
	new("SosObservationOffering", id = id, name = name,
			time = time, procedure = procedure,
			observedProperty = observedProperty,
			featureOfInterest = featureOfInterest,
			responseFormat = responseFormat,
			intendedApplication = intendedApplication,
			resultModel = resultModel, responseMode = responseMode,
			boundedBy = boundedBy)
}

SosContents <- function(observationOfferings) {
	new("SosContents", observationOfferings = observationOfferings)
}

SosEventTime <- function(temporalOps) {
	new("SosEventTime", temporalOps = temporalOps)
}

SosEventTimeLatest <- function() {
	new("SosEventTimeLatest")
}

SosFeatureOfInterest <- function(objectIDs = list(NA), spatialOps = NULL) {
	new("SosFeatureOfInterest", objectIDs = objectIDs, spatialOps = spatialOps)
}

#
#
#
SosDescribeSensor <- function(
		service,
		version,
		procedure,
		outputFormat) {
	new("SosDescribeSensor",
			request = sosDescribeSensorName,
			service = service,
			version = version,
			procedure = procedure,
			outputFormat = outputFormat)
}

#
#
#
SosGetObservation <- function(
		service,
		version,
		offering, 
		observedProperty,
		responseFormat, 
		srsName = as.character(NA),
		eventTime = list(NA), 
		procedure = as.character(NA),
		featureOfInterest = NULL, 
		result = NULL,
		resultModel = as.character(NA),
		responseMode = as.character(NA),
		BBOX = as.character(NA)) {
	new("SosGetObservation",
			request = sosGetObservationName,
			service = service,
			version = version,
			offering = offering, 
			observedProperty = observedProperty,
			responseFormat = responseFormat,
			srsName = srsName,
			eventTime = eventTime,
			procedure = procedure, 
			featureOfInterest = featureOfInterest,
			result = result,
			resultModel = resultModel,
			responseMode = responseMode,
			BBOX = BBOX)
}

#
#
#
SosGetObservationById <- function(
		service,
		version,
		observationId,
		responseFormat, 
		srsName = as.character(NA),
		resultModel = as.character(NA),
		responseMode = as.character(NA)) {
	new("SosGetObservationById",
			request = sosGetObservationByIdName,
			service = service,
			version = version,
			observationId = observationId,
			responseFormat = responseFormat,
			srsName = srsName,
			resultModel = resultModel,
			responseMode = responseMode)
}


################################################################################
# main request method
#
.sosRequest_1.0.0 <- function(sos, request, verbose = FALSE, inspect = FALSE) {
	# check the request for consistency with service description
	.checkResult <- checkRequest(service = sos, operation = request,
			verbose = verbose)
	if(!.checkResult) {
		warning("Check returned FALSE! Turn on verbose option for possible details.",
				immediate. = TRUE)
	}
	
	.response = ""
	
	# get encoding function for the respective method
	.encodingFunction <- sos@encoders[[sos@method]]
	if(verbose) {
		.f <- functionBody(.encodingFunction)
		cat("[.sosRequest_1.0.0] Encoding Function (beginning of function body): ",
				substring(text = .f, first = 0, last = 60), " ... [",
				max((length(.f) - 60), 0), " more chrs].\n")
	}
	
	# encode!
	.encodedRequest = .encodingFunction(obj = request, sos = sos,
			verbose = verbose)
	
	if(sos@method == .sosConnectionMethodGet) {
		.dcp <- sosGetDCP(sos, sosName(request), "Get") #sos@url
		if(is.null(.dcp) || is.na(.dcp)) {
			.dcp <- sos@url
			if(verbose) cat("[.sosRequest_1.0.0] Could not get DCP from operation description. This is OK for first GetCapabilities request.\n")
		}
			
		if(isTRUE(grep(pattern = "[\\?]", x = .dcp) > 0)) {
			warning("Given url already contains a '?', appending arguments!")
			.url = paste(.dcp, .encodedRequest, sep = "&")
		}
		else .url = paste(.dcp, .encodedRequest, sep = "?")
		
		if(verbose || inspect) {
			cat("[.sosRequest_1.0.0] GET!\n[.sosRequest_1.0.0] REQUEST: ", .url,
					"\n")
		}
		
		.response = getURL(url = .url, .opts = sos@curlOptions,
				curl = sos@curlHandle,
				.encoding = sosDefaultCharacterEncoding)
	}
	else if(sos@method == .sosConnectionMethodPost) {
		if(verbose || inspect) {
			cat("[.sosRequest_1.0.0] POST!\n[.sosRequest_1.0.0] REQUEST:")
			print(.encodedRequest)
		}
		
		.dcp <- sosGetDCP(sos, sosName(request), "Post") #sos@url
		if(is.null(.dcp) || is.na(.dcp)) {
			.dcp <- sos@url
			if(verbose) cat("[.sosRequest_1.0.0] Could not get DCP from operation description.  This is OK for first GetCapabilities request.\n")
		}
		
		# using 'POST' for application/x-www-form-urlencoded content
		.response <- postForm(uri = .dcp,
				request = toString(.encodedRequest),
				style = "POST", .opts = sos@curlOptions,
				curl = sos@curlHandle,
				.encoding = sosDefaultCharacterEncoding)
	}
	else if(sos@method == .sosConnectionMethodSOAP) {
		if(verbose || inspect) {
			print("[.sosRequest_1.0.0] SOAP! REQUEST:\n")
			print(.encodedRequest)
		}
		
		# TODO add SOAP request method
	}
	else {
		stop(paste("Unsupported method, has to be one of",
						SosSupportedConnectionMethods()))
	}
	
	if(verbose) {
		cat("[.sosRequest_1.0.0] response:\n")
		print(.response)
		if(is.raw(.response)) cat("raw as char: ", rawToChar(.response), "\n")
	}
	
	if(length(.response) > 0 & 
			regexpr("(<html>|<HTML>|<!DOCTYPE HTML)", .response) > 0) {
		if(verbose) cat("[.sosRequest_1.0.0] Got HTML, probably an error.\n")
		
		# might still be KML with embedded HTML!
		if(regexpr("(http://www.opengis.net/kml/)", .response) > 0) {
			if(verbose) cat("[.sosRequest_1.0.0] Got KML! Can continue...\n")
		}
		else stop(paste("[sos4R] ERROR: Got HTML response!:\n", .response,
							"\n\n"))
	}
	
	return(.response)
}

setMethod(f = "sosRequest",
		signature = signature(sos = "SOS_1.0.0", request = "OwsServiceOperation",
				verbose = "logical", inspect = "logical"),
		def = function(sos, request, verbose, inspect) {
			.sosRequest_1.0.0(sos = sos, request = request, verbose = verbose,
					inspect = inspect)
		}
)


################################################################################
# functions for SOS operations

#
#
#
.getCapabilities_1.0.0 <- function(sos, verbose, inspect, sections,
		acceptFormats, updateSequence, owsVersion,	acceptLanguages) {
	if (verbose) {
		cat("[.getCapabilities_1.0.0] of", sosUrl(sos), "\n")
	}
	
	.gc <- OwsGetCapabilities(service = sosService,
			acceptVersions = c(sosVersion(sos)), sections = sections,
			acceptFormats = acceptFormats, updateSequence = updateSequence,
			owsVersion = owsVersion, acceptLanguages = acceptLanguages)
	if(verbose) cat("[.getCapabilities_1.0.0] REQUEST:\n", toString(.gc), "\n")
	
	.responseString = sosRequest(sos = sos, request = .gc,
			verbose = verbose, inspect = inspect)
	if(verbose){
		cat("[.getCapabilities_1.0.0] RESPONSE:\n", .responseString , "\n")
	}
	
	.response <- xmlParseDoc(file = .responseString, asText = TRUE)
	if(verbose || inspect) {
		cat("[.getCapabilities_1.0.0] RESPONSE DOC:\n")
		print(.response)
	}
	
	if(.isExceptionReport(.response)) {
		return(.handleExceptionReport(sos, .response))
	}
	else {
		.parsingFunction <- sosParsers(sos)[[sosGetCapabilitiesName]]
		.caps <- .parsingFunction(obj = .response, sos = sos)
		if (verbose) {
			cat("[.getCapabilities_1.0.0] DONE WITH PARSING!\n")
		} 
		return(.caps)
	}
}
setMethod(f = "getCapabilities", signature = signature(sos = "SOS_1.0.0"),
		def = function(sos, verbose, inspect, sections, acceptFormats,
				updateSequence, owsVersion,	acceptLanguages) {
			return(.getCapabilities_1.0.0(sos = sos, verbose = verbose,
							inspect = inspect, sections = sections,
							acceptFormats = acceptFormats,
							updateSequence = updateSequence,
							owsVersion = owsVersion,
							acceptLanguages = acceptLanguages))
		}
)

#
#
#
.describeSensor_1.0.0 <- function(sos, procedure, outputFormat, verbose,
		inspect, saveOriginal) {
	if(verbose) cat("[.describeSensor_1.0.0] ", procedure, "@", sos@url, "\n")
	
	# check if multiple sensors
	if(length(procedure) > 1) {
		if(verbose) cat("[.describeSensor_1.0.0] multiple sensors: ", procedure,
					"\n")
		
		.descriptions <- list()
		for (p in procedure) {
			.description <- .describeSensor_1.0.0(sos = sos, procedure = p,
					outputFormat = outputFormat, verbose = verbose,
					inspect = inspect, saveOriginal = saveOriginal)
			.descriptions <- c(.descriptions, .description)
		}
		
		return(.descriptions)
	}
	
	.ds <- SosDescribeSensor(service = sosService, version = sos@version,
			procedure = procedure, outputFormat = outputFormat)
	if(verbose)
		cat("[.describeSensor_1.0.0] REQUEST:\n", toString(.ds), "\n")
	
	
	.responseString = sosRequest(sos = sos, request = .ds,
			verbose = verbose, inspect = inspect)
	if(verbose || inspect){
		cat("[.describeSensor_1.0.0] RESPONSE:\n", .responseString , "\n")
	}
	
	.response <- xmlParseDoc(.responseString, asText = TRUE)
	if(verbose || inspect) {
		cat("[.describeSensor_1.0.0] RESPONSE DOC:\n")
		print(.response)
	}
	
	.filename <- NULL
	if(!is.null(saveOriginal)) {
		if(is.character(saveOriginal)) {
			.filename <- paste(saveOriginal, ".xml", sep = "")
			if(verbose) cat("Using saveOriginal parameter for file name:",
						.filename, "\n")
		} 
		else if(is.logical(saveOriginal)) {
			if(saveOriginal) .filename <- paste(.cleanupFileName(procedure),
						".xml", sep = "")
			if(verbose) cat("Generating file name:", .filename, "\n")
		}
		
		if(verbose) {
			cat("[.describeSensor_1.0.0] Saving original document...",
					.filename, "in", getwd(), "\n")
		}
		
		# TODO alternatively one could use tempfile() instead of implicit getwd()
		saveXML(.response, file = .filename)
		
		cat("[sos4R] Original document saved:", .filename, "\n")
	}
	
	if(.isExceptionReport(.response)) {
		return(.handleExceptionReport(sos, .response))
	}
	else {
		.parsingFunction <- sosParsers(sos)[[sosDescribeSensorName]]
		.sml <- .parsingFunction(obj = .response, sos = sos, verbose = verbose)
		
		if(!is.null(.filename)) {
			.oldAttrs <- attributes(.sml)
			.newAttrs <- list(.filename)
			names(.newAttrs) <- list(sosAttributeFileName)
			if(verbose) cat("[.describeSensor_1.0.0] Appending new attributes",
						toString(.newAttrs), "(names",
						toString(names(.newAttrs)), ")\n")
			
			attributes(.sml) <- c(.oldAttrs, .newAttrs)
		}
		
		return(.sml)
	}
}
setMethod(f = "describeSensor",
		signature = signature(sos = "SOS_1.0.0", procedure  = "character"), 
		def = function(sos, procedure, outputFormat, verbose, inspect,
				saveOriginal) {
			.result <- .describeSensor_1.0.0(sos = sos, procedure = procedure,
					outputFormat = outputFormat, verbose = verbose,
					inspect = inspect, saveOriginal = saveOriginal)
			return(.result)
		}
)


#
# 
#
setMethod(f = "getObservationById",
		signature = signature(sos = "SOS_1.0.0", observationId = "character"), 
		def = function(sos, observationId, responseFormat, srsName,
				resultModel, responseMode, verbose, inspect, saveOriginal) {
			return(.getObservationById_1.0.0(sos = sos,
							observationId = observationId,
							responseFormat = responseFormat, srsName = srsName,
							resultModel = resultModel,
							responseMode = responseMode, verbose = verbose,
							inspect = inspect, saveOriginal = saveOriginal))
		}
)

.getObservationById_1.0.0 <- function(sos, observationId, responseFormat, srsName,
		resultModel, responseMode, verbose, inspect, saveOriginal) {
	if(verbose) {
		cat("[.getObservationById_1.0.0] ID", observationId, "\n")
	}
	
	.filename <- NULL
	if(is.character(saveOriginal)) {
		.filename <- saveOriginal
		if(verbose) cat("[.getObservationById_1.0.0] Using saveOriginal parameter for file name:",
					.filename, "\n")
	} 
	else if(is.logical(saveOriginal)) {
		if(saveOriginal) .filename <- paste(observationId, 
					format(Sys.time(), sosDefaultFilenameTimeFormat),
					sep = "_")
		if(verbose) cat("[.getObservationById_1.0.0] Generating file name:", .filename, "\n")
	}
	
	.go <- SosGetObservationById(service = sosService,
			version = sos@version, observationId = observationId,
			responseFormat =  responseFormat, srsName = srsName,
			resultModel = resultModel, responseMode = responseMode)
	
	if(verbose)
		cat("[.getObservationById_1.0.0] REQUEST:\n", toString(.go), "\n")
	
	.responseString = sosRequest(sos = sos, request = .go,
			verbose = verbose, inspect = inspect)
	if(verbose || inspect){
		cat("[.getObservationById_1.0.0] RESPONSE:\n", .responseString , "\n")
	}
	
	.response <- xmlParseDoc(.responseString, asText = TRUE)
	if(verbose || inspect) {
		cat("[.getObservationById_1.0.0] RESPONSE DOC:\n")
		print(.response)
	}
	
	if(!is.null(.filename)) {
		.filename <- paste(.filename, ".xml", sep = "")
		saveXML(.response, file = .filename)
		cat("[sos4R] Original document saved:", .filename, "\n")
	}
	
	if(.isExceptionReport(.response)) {
		return(.handleExceptionReport(sos, .response))
	}
	else {
		.parsingFunction <- sosParsers(sos)[[sosGetObservationByIdName]]
		.obs <- .parsingFunction(obj = .response, sos = sos,
				verbose = verbose)
		
		# remove list if only one element
		if(is.list(.obs) && length(.obs) == 1)
			.obs <- .obs[[1]]
		
		if(verbose) {
			cat("[.getObservationById_1.0.0] PARSED RESPONSE:\n")
			print(.obs)
		}
		
		if(!is.null(.filename)) {
			.oldAttrs <- attributes(.obs)
			.newAttrs <- list(.filename)
			names(.newAttrs) <- list(sosAttributeFileName)
			if(verbose) cat("[.getObservationById_1.0.0] Appending new attributes",
						toString(.newAttrs), "(names",
						toString(names(.newAttrs)), ")\n")
			
			attributes(.obs) <- c(.oldAttrs, .newAttrs)
		}
		
		return(.obs)
	}
	
	if(verbose) {
		cat("[.getObservationById_1.0.0] returning raw response string.\n")
	}
	
	return(.responseString)
}


#
#
#
.getObservation_1.0.0 <- function(sos, offeringId, observedProperty,
		responseFormat, srsName, eventTime,	procedure, featureOfInterest,
		result, resultModel, responseMode, BBOX, latest, verbose, inspect,
		saveOriginal) {
	
	.filename <- NULL
	if(is.character(saveOriginal)) {
		.filename <- saveOriginal
		if(verbose) cat("[.getObservation_1.0.0] Using saveOriginal parameter for file name:",
					.filename, "\n")
	} 
	else if(is.logical(saveOriginal)) {
		if(saveOriginal) .filename <- paste(.cleanupFileName(offeringId), 
					format(Sys.time(), sosDefaultFilenameTimeFormat), sep = "_")
		if(verbose) cat("[.getObservation_1.0.0] Generating file name:",
					.filename, "\n")
	}
	
	if(verbose)
		cat("[.getObservation_1.0.0] to ", sos@url, " with offering ",
				offeringId, "\n")
	
	if(latest) .eventTime <- list(.createLatestEventTime(verbose))
	else .eventTime <- eventTime
	
	if(latest && !is.na(eventTime))
		warning("'Latest' is set to TRUE > given eventTime is ignored!")
	
	.go <- SosGetObservation(service = sosService, version = sos@version, 
			offering = offeringId, observedProperty = observedProperty,
			responseFormat =  responseFormat, srsName = srsName,
			eventTime = .eventTime, procedure = procedure,
			featureOfInterest = featureOfInterest, result = result,
			resultModel = resultModel, responseMode = responseMode,
			BBOX = BBOX)
	
	if(verbose)
		cat("[.getObservation_1.0.0] REQUEST:\n", toString(.go), "\n")
	
	.responseString = sosRequest(sos = sos, request = .go,
			verbose = verbose, inspect = inspect)
	
	cat("[sos4R] Received response (size:", object.size(.responseString),
			"bytes), parsing ...\n")
	
	# responseFormat starts with text/xml OR the response string is XML content,
	# for example an exeption (which is xml even if request wants something else
	.contentType <- NA_character_
	.contentType <- attributes(.responseString)[["Content-Type"]]
	
	if(verbose) cat("[.getObservation_1.0.0] Content-Type:", .contentType, "\n")
	
	if(isXMLString(.responseString)) {
		if(verbose) cat("[.getObservation_1.0.0] Got XML string as response",
					"(based on isXMLString()).\n")
		
		.hasSubtype <- FALSE
		.contentSubtype <- NA
		if(length(.contentType) < 1) {
			if(verbose) cat("[.getObservation_1.0.0] No content type!",
						"Falling back to '", mimeTypeXML, "'\n")
			.contentType <- mimeTypeXML
		}
		else if(length(.contentType) > 1) {
			# check if subtype is present or take just the first
			.subtypeIdx <- which(names(.contentType) == "subtype")
			if(.subtypeIdx > 0) {
				.hasSubtype <- TRUE
				.contentSubtype <- .contentType[[.subtypeIdx]]
				if(verbose) cat("[.getObservation_1.0.0] Found mime subtype: ",
							toString(.contentSubtype), "'\n")
			}
			else if(verbose) cat(
						"[.getObservation_1.0.0] More than one content type, ",
								"no subtype detected : '",
								toString(.contentType),
								"'\n\tUsing the first one: '",
								.contentType[[1]], "'\n")
			.contentType <- .contentType[[1]]
		}

		.response <- xmlParseDoc(.responseString, asText = TRUE)
		if(verbose || inspect) {
			cat("[.getObservation_1.0.0] RESPONSE DOC:\n")
			print(.response)
		}
		# select the parser and file ending based on the mime type FIRST
		.fileEnding <- ".xml"
		if(.contentType == mimeTypeXML) {
			if(.hasSubtype && .contentSubtype == mimeSubtypeOM) {
				if(verbose)
					cat("[.getObservation_1.0.0] Got OM according to mime type.\n")
				.parserName <- mimeTypeOM
			}
			else {
				if(verbose)
					cat("[.getObservation_1.0.0] Got pure XML according to mime type. ",
							"Trying to parse with default parser, see SosParsingFunctions().\n")
				.parserName <- mimeTypeXML
			}
		}
		else if (.contentType == mimeTypeKML) {
			if(verbose) cat("[.getObservation_1.0.0] Got KML according to mime type.\n")
			
			.fileEnding <- ".kml"
			.parserName <- mimeTypeKML
		}
		else {
			# fall back, or more of a default: the function name
			.parserName <- sosGetObservationName
		}
		
		if(!is.null(.filename)) {
			.filename <- paste(.filename, .fileEnding, sep = "")
			saveXML(.response, file = .filename)
			
			if(verbose) {
				cat("[.getObservation_1.0.0] Saved original document:",
						.filename, "\n")
			}
		}
		
		if(.isExceptionReport(.response)) {
			return(.handleExceptionReport(sos, .response))
		}
		
		if( !is.na(responseFormat) && 
				isTRUE(grep(pattern = "text/xml", x = responseFormat) != 1)) {
			warning("Got XML string, but request did not require text/xml (or subtype).")
		}
		
		.parsingFunction <- sosParsers(sos)[[.parserName]]
		
		if(verbose) {
			cat("[.getObservation_1.0.0] Parsing with function ")
			print(.parsingFunction)
		}
		
		.obs <- .parsingFunction(obj = .response, sos = sos,
				verbose = verbose)
		
		# calculate result length vector
		if(inherits(.obs, "OmObservationCollection")) {
			.resultLength <- sapply(sosResult(.obs, bind = FALSE,
							coordinates = FALSE), nrow)
			if(length(.resultLength) == 0) # nothing
				.resultLength = 0
		}
		else .resultLength <- NA
		
		if(verbose) {
			cat("[.getObservation_1.0.0] PARSED RESPONSE:",
					class(.obs), "\n")
			cat("[.getObservation_1.0.0] Result length(s): ",
					toString(.resultLength), "\n")
		}
		
		if(is.list(.obs) && any(sapply(.obs, is.null))) {
			.countInfo <- paste("NO DATA, turn on 'verbose' for more information.")
		}
		else {
			.countInfo <- paste(sum(.resultLength), "result values [",
					toString(.resultLength), "].")
		}

		.msg <- paste("[sos4R] Finished getObservation to", sos@url,
				"\n\t--> received", length(.obs), "observation(s) having",
				.countInfo , "\n")
		if(!is.null(.filename)) {
			.msg <- paste(.msg,
					"[sos4R] Original document saved:", .filename, "\n")
			
			.oldAttrs <- attributes(.obs)
			.newAttrs <- list(.filename)
			names(.newAttrs) <- list(sosAttributeFileName)
			if(verbose) cat("[.getObservationById_1.0.0] Appending new attributes",
						toString(.newAttrs), "(names",
						toString(names(.newAttrs)), ")\n")
			
			attributes(.obs) <- c(.oldAttrs, .newAttrs)
		}
		cat(.msg)
		
		# RETURN ###
		return(.obs)
	}
	else { # response is NOT an XML string:
		if(verbose)
			cat("[.getObservation_1.0.0] Did NOT get XML string as response, trying to parse with",
					responseFormat, "\n")
		
		if(mimeTypeCSV == responseFormat) {
			if(verbose || inspect) {
				cat("[.getObservation_1.0.0] CSV RESPONSE:\n")
				print(.responseString)
			}
		
			.parsingFunction <- sosParsers(sos)[[mimeTypeCSV]]
			.csv <- .parsingFunction(obj = .responseString, verbose = verbose)
		
			if(!is.null(.filename)) {
				.filename <- paste(file = .filename, ".csv", sep = "")
				write.csv(.csv, .filename)
			}
			
			.msg <- paste("[sos4R] Finished getObservation to", sos@url, "\n\t",
					"--> received observations with dimensions", 
					toString(dim(.csv)), "\n")
			if(!is.null(.filename)) {
				.msg <- paste(.msg,
						"[sos4R] Original document saved:", .filename, "\n")
				
				.oldAttrs <- attributes(.csv)
				.newAttrs <- list(.filename)
				names(.newAttrs) <- list(sosAttributeFileName)
				if(verbose) cat("[.getObservation_1.0.0] Appending new attributes",
							toString(.newAttrs), "(names",
							toString(names(.newAttrs)), ")\n")
				
				attributes(.csv) <- c(.oldAttrs, .newAttrs)
			}
			cat(.msg)
		
			# RETURN ###
			return(.csv)
		} # grep(pattern = mimeTypeCSV...
		
		# Add other non-XML encodings here.
	} # else

	# not xml nor csv
	if(verbose || inspect) {
		cat("[.getObservation_1.0.0] UNKNOWN RESPONSE FORMAT:\n")
		cat(.responseString, "\n")
		cat("[.getObservation_1.0.0] Content-Type: ", .contentType)
		warning("Unknown response format!")
	}
	
	if(!is.null(.filename)) {
		save(.responseString, file = .filename)
		cat("[sos4R] Saved original document:", .filename)
	}
	else warning("File name is NULL, could not save document!")
	
	# RETURN ##############
	return(.responseString)
}

#
#
#
setMethod(f = "getObservation",
		signature = signature(sos = "SOS_1.0.0",
				offering = "SosObservationOffering"),
		def = function(sos, offering, observedProperty, responseFormat, srsName,
				eventTime,	procedure, featureOfInterest, result, resultModel,
				responseMode, BBOX, latest, verbose, inspect, saveOriginal) {
			.offeringId <- offering@id
			if(verbose)	cat("[getObservation] Requesting offering", .offeringId,
						"by SosObservationOffering.\n")
			
			return(.getObservation_1.0.0(sos = sos, offeringId = .offeringId,
							observedProperty = observedProperty,
							responseFormat = responseFormat,
							srsName = srsName, eventTime = eventTime,
							procedure = procedure,
							featureOfInterest = featureOfInterest,
							result = result, resultModel = resultModel,
							responseMode = responseMode, BBOX = BBOX,
							latest = latest, verbose = verbose,
							inspect = inspect, saveOriginal = saveOriginal))
		}
)

#
#
#
setMethod(f = "getObservation",
		signature = signature(sos = "SOS_1.0.0",
				offering = "character"),
		def = function(sos, offering, observedProperty = list(), responseFormat,
				srsName, eventTime,	procedure, featureOfInterest, result,
				resultModel, responseMode, BBOX, latest, verbose, inspect,
				saveOriginal) {
			if(verbose)	cat("[getObservation] Requesting offering", offering,
						"by name.\n")
			
			.off <- sosOfferings(sos)[[offering]]
			
			if(length(observedProperty) == 0) {
				.obsProps <- sosObservedProperties(.off)
				if(verbose) cat("[getObservation] Got observation(s) from offering because none given:",
							toString(.obsProps), "\n")
			}
			else {
				.obsProps <- observedProperty
			}
			
			return(.getObservation_1.0.0(sos = sos, offeringId = offering,
							observedProperty = .obsProps,
							responseFormat = responseFormat,
							srsName = srsName, eventTime = eventTime,
							procedure = procedure,
							featureOfInterest = featureOfInterest,
							result = result, resultModel = resultModel,
							responseMode = responseMode, BBOX = BBOX,
							latest = latest, verbose = verbose,
							inspect = inspect, saveOriginal = saveOriginal))
		}
)

#
# see: http://www.oostethys.org/best-practices/best-practices-get
#
setMethod("encodeRequestKVP", "SosDescribeSensor", 
		function(obj, sos, verbose = FALSE) {
			
			if(obj@version == "1.0.0") {
				return(.sosEncodeRequestKVPDescribeSensor_1.0.0(obj = obj,
								sos = sos, verbose = verbose))
			}
			else {
				stop("Version not supported!")
			}
		}
)
.sosEncodeRequestKVPDescribeSensor_1.0.0 <- function(obj, sos,
		verbose = FALSE) {
	# mandatory:
	.service <- paste("service",
			.kvpEscapeSpecialCharacters(x = obj@service), sep = "=")
	.request <- paste("request" , sosDescribeSensorName, sep = "=")
	.version <- paste("version", 
			.kvpEscapeSpecialCharacters(x = obj@version), sep = "=")
	.procedure <- paste("procedure",
			.kvpEscapeSpecialCharacters(x = obj@procedure), sep = "=")
	.format <- paste(
			"outputFormat",
			.kvpEscapeSpecialCharacters(x = gsub(obj@outputFormat,
							pattern = "&quot;",
							replacement = '"')),
			sep = "=")
	
	.kvpString <- paste(.service, .request, .version, .procedure,
			.format, sep = "&")
	
	if(verbose)
		cat(.kvpString)
	
	return(.kvpString)
}
setMethod("encodeRequestKVP", "SosGetObservation", 
		function(obj, sos, verbose = FALSE) {
			if(obj@version == "1.0.0") {
				return(.sosEncodeRequestKVPGetObservation_1.0.0(obj, sos,
								verbose))		
			}
			else {
				stop("Version not supported!")
			}
		}
)
.sosEncodeRequestKVPGetObservation_1.0.0 <- function(obj, sos,
		verbose = FALSE) {
	if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] encoding",
				toString(obj))
	
	# required:
	.request <- paste("request" , sosGetObservationName, sep = "=")
	.service <- paste("service",
			.kvpEscapeSpecialCharacters(x = obj@service), sep = "=")
	.version <- paste("version",
			.kvpEscapeSpecialCharacters(x = obj@version), sep = "=")
	.offering <- paste("offering",
			.kvpEscapeSpecialCharacters(x = obj@offering), sep = "=")
	.observedProperty <- .kvpKeyAndValues("observedProperty", 
			obj@observedProperty)
	
	.mandatory <- paste(.service, .request, .version, .offering,
			.observedProperty, sep = "&")
	
	if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0]",
				"mandatory elements: ", .mandatory)
	
	# optional:
	.optionals = ""
	# is optional for GET
	if( !is.na(obj@responseFormat)) {
		if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] Adding response format ",
					obj@responseFormat, "\n")
		.responseFormat <- paste(
				"responseFormat", 
				.kvpEscapeSpecialCharacters(x = gsub(obj@responseFormat,
								pattern = "&quot;",
								replacement = '"')),
				sep = "=")
		.optionals <- paste(.optionals, .responseFormat, sep = "&")
	}
	
	if( !is.na(obj@srsName)) {
		if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] Adding SRS name ",
					obj@srsName, "\n")
		.optionals <- paste(.optionals, paste("srsName", 
						.kvpEscapeSpecialCharacters(x = obj@srsName),
						sep = "="),
				sep = "&")
	}
	
	if( !is.na(obj@eventTime)) {
		if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] Adding event time",
					toString(obj@eventTime), "\n")
		if(length(obj@eventTime) > 1)
			warning("Only first event time in the list is used for KVP!")
		
		.timeString <- encodeKVP(obj = obj@eventTime[[1]],
				sos = sos, verbose = verbose)
		
		# if the eventTime is a latest request, it returns NA, the GET binding
		# says for the latest observation eventTime is omitted
		if(!is.na(.timeString)) {
			.optionals <- paste(.optionals, paste("eventTime", 
							.kvpEscapeSpecialCharacters(x = .timeString), 
							sep = "="), 
					sep = "&")
		}
		else {
			if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] ", 
						"encodeKVP returned NA for eventTime, omitting",
						"parameter for request for latest observation.")
		}
	}
	
	if( !any(sapply(obj@procedure, "is.na"), na.rm = TRUE)) {
		if(verbose) cat("Adding procedures ", obj@procedure, "\n")
		.optionals <- paste(.optionals, .kvpKeyAndValues("procedure",
						obj@procedure), sep = "&")
	}
	
	if( !is.null(obj@featureOfInterest)) {
#		print(obj@featureOfInterest)
#		.optionals <- paste(.optionals, .kvpKeyAndValues("featureOfInterest",
#						obj@featureOfInterest), sep = "&")
		warning("'featureOfInterest' is not supported for 'GET' - parameter is discarded, use another method to include it!")
	}
	
	if( !is.null(obj@result)) {
		warning("'result' is not supported for 'GET' - parameter is discarded, use another method to include it!")
	}
	
	if( !is.na(obj@resultModel)) {
		if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] Adding result model ",
					obj@resultModel, "\n")
		.optionals <- paste(.optionals, paste("resultModel",
						.kvpEscapeSpecialCharacters(x = obj@resultModel),
						sep = "="),
				sep = "&")
	}
	
	if( !is.na(obj@responseMode)) {
		if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] Adding response mode ",
					obj@responseMode, "\n")
		.optionals <- paste(.optionals, paste("responseMode",
						.kvpEscapeSpecialCharacters(x = obj@responseMode),
						sep = "="),
				sep = "&")
	}
	
	if( !is.na(obj@BBOX)) {
		if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0] Adding BBOX ",
					obj@BBOX, "\n")
		.optionals <- paste(.optionals, paste("BBOX", 
						.kvpEscapeSpecialCharacters(x = obj@BBOX), sep = "="),
				sep = "&")
	}
	
	if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0]",
				"optional elements: ", .optionals)
	
	.kvpString <- paste(.mandatory, .optionals, sep = "")
	
	if(verbose) cat("[.sosEncodeRequestKVPGetObservation_1.0.0]",
				"Finished KVP string creation:\n", .kvpString, "\n")
	
	return(.kvpString)
}

setMethod("encodeRequestKVP", "SosGetObservationById", 
		function(obj, sos, verbose = TRUE) {
			stop("KVP encoding of operation 'GetObservationById' not supported!")
		}
)

#
# encode as XML
#
setMethod("encodeRequestXML", "SosGetObservation", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeRequestXML]", class(obj), "\n")
			}
			
			if(obj@version == "1.0.0") {
				return(.sosEncodeRequestXMLGetObservation_1.0.0(obj = obj,
								sos = sos, verbose = verbose))		
			}
			else {
				stop("Version not supported!")
			}
		}
)
.sosEncodeRequestXMLGetObservation_1.0.0 <- function(obj, sos,
		verbose = FALSE) {
	.xmlDoc <- xmlNode(name = sosGetObservationName,
			namespace = sosNamespacePrefix,
			namespaceDefinitions = c(.sosNamespaceDefinitionsForAll,
					.sosNamespaceDefinitionsGetObs),
			attrs=c(.xsiSchemaLocationAttribute, service = obj@service,
					version = obj@version))
	
	# required and optional are mixed - schema requires a particular order:
	.offering <- xmlNode(name = "offering", namespace = sosNamespacePrefix,
			obj@offering)
	.xmlDoc <- addChildren(node = .xmlDoc, .offering)
	
	if(!any(is.na(obj@eventTime))) {
		.eventTimeList <- lapply(obj@eventTime, encodeXML, sos = sos,
				verbose = verbose)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = .eventTimeList,
				append = TRUE)
	}
	
	if( !any(sapply(obj@procedure, "is.na"), na.rm = TRUE)) {
		.procedureList <- lapply(obj@procedure, "xmlNode",
				name="procedure", namespace = sosNamespacePrefix)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = .procedureList,
				append = TRUE)
	}
	
	.observedProperties <- lapply(obj@observedProperty, "xmlNode",
			name="observedProperty", namespace = sosNamespacePrefix)
	.xmlDoc <- addChildren(node = .xmlDoc, kids = .observedProperties,
			append = TRUE)
	
	if( !is.null(obj@featureOfInterest)) {
		.foi <- encodeXML(obj = obj@featureOfInterest, sos = sos,
				verbose = verbose)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.foi),
				append = TRUE)
	}
	
	if( !is.null(obj@result)) {
		if(is.character(obj@result)) {
			.result <- encodeXML(obj = obj@result, sos = sos,
					addNamespaces = TRUE, verbose = verbose)
		}
		else {
			.result <- encodeXML(obj = obj@result, sos = sos, verbose = verbose)
		}
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.result),
				append = TRUE)
	}
	
	if( !is.na(obj@responseFormat)) {
		.rF <- gsub(obj@responseFormat, pattern = "&quot;", replacement = "\"")
		
		.responseFormat <- xmlNode(name = "responseFormat",
				namespace = sosNamespacePrefix, value = .rF)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.responseFormat),
				append = TRUE)
	}
	
	if( !is.na(obj@resultModel)) {
		.resultModel <- xmlNode(name = "resultModel",
				namespace = sosNamespacePrefix,
				obj@resultModel)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.resultModel),
				append = TRUE)
	}
	
	if( !is.na(obj@responseMode)) {
		.responseMode <- xmlNode(name = "responseMode",
				namespace = sosNamespacePrefix,
				obj@responseMode)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.responseMode),
				append = TRUE)
	}
	
	if( !is.na(obj@srsName)) {
		.xmlDoc <- addAttributes(.xmlDoc, srsName = obj@srsName, append = TRUE)
	}
	
	if( !is.na(obj@BBOX)) {
		warning("GetObservation contains BBOX, but that is not supported for 'POST' (and not at all in the SOS Specification...) - use featureOfInterest instead!")
	}
	
	return(.xmlDoc)
}

setMethod("encodeRequestXML", "SosGetObservationById", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeRequestXML]", class(obj), "\n")
			}
			
			if(obj@version == "1.0.0") {
				return(.sosEncodeRequestXMLGetObservationById_1.0.0(obj = obj,
								sos = sos))		
			}
			else {
				stop("Version not supported!")
			}
		}
)
.sosEncodeRequestXMLGetObservationById_1.0.0 <- function(obj, sos) {
	.xmlDoc <- xmlNode(name = "GetObservationById",
			namespace = sosNamespacePrefix,
			namespaceDefinitions = c(.sosNamespaceDefinitionsForAll,
					.sosNamespaceDefinitionsGetObs),
			attrs=c(.xsiSchemaLocationAttribute,
					service = obj@service, version = obj@version))
	
	.obsId <- xmlNode(name = "ObservationId", namespace = sosNamespacePrefix,
			obj@observationId)
	.xmlDoc <- addChildren(node = .xmlDoc, .obsId)
	
	.rF <- gsub(obj@responseFormat, pattern = "&quot;", replacement = "\"")
	.responseFormat <- xmlNode(name = "responseFormat",
			namespace =  sosNamespacePrefix, .rF)
	.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.responseFormat),
			append = TRUE)
	
	if( !is.na(obj@resultModel)) {
		.resultModel <- xmlNode(name = "resultModel",
				namespace =  sosNamespacePrefix,
				obj@resultModel)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.resultModel),
				append = TRUE)
	}
	
	if( !is.na(obj@responseMode)) {
		.responseMode <- xmlNode(name = "responseMode",
				namespace =  sosNamespacePrefix,
				obj@responseMode)
		.xmlDoc <- addChildren(node = .xmlDoc, kids = list(.responseMode),
				append = TRUE)
	}
	
	if( !is.na(obj@srsName)) {
		.xmlDoc <- addAttributes(.xmlDoc, srsName = obj@srsName, append = TRUE)
	}
	
	return(.xmlDoc)
}


#
# encode as XML
#
setMethod("encodeRequestXML", "SosDescribeSensor", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeRequestXML]", class(obj), "\n")
			}
			
			if(obj@version == "1.0.0") {
				return(.sosEncodeRequestXMLDescribeSensor_1.0.0(obj = obj))
			}
			else {
				stop("Version not supported!")
			}
		}
)
.sosEncodeRequestXMLDescribeSensor_1.0.0 <- function(obj) {
	xmlDoc <- xmlNode(name = sosDescribeSensorName,
			namespace = sosNamespacePrefix,
			namespaceDefinitions = .sosNamespaceDefinitionsForAll,
			attrs=c(.xsiSchemaLocationAttribute,
					service = obj@service,
					outputFormat = obj@outputFormat,
					version = obj@version))
	
	procedure <- xmlNode(name = "procedure", namespace = sosNamespacePrefix,
			obj@procedure)
	xmlDoc$children[[1]] <- procedure
	
	return(xmlDoc)
}

#
# encode for SOAP
#
setMethod("encodeRequestSOAP", "SosDescribeSensor", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("ENCODE SOAP ", class(obj), "\n")
			}
			
			if(obj@version == "1.0.0") {
				return(.sosEncodeRequestXMLDescribeSensor_1.0.0(obj))
			}
			else {
				stop("Version not supported!")
			}
		}
)
setMethod("encodeRequestSOAP", "SosGetObservation", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("ENCODE SOAP ", class(obj), "\n")
			}
			stop("Method not implemented yet!")
		}
)
setMethod("encodeRequestSOAP", "SosGetObservationById", 
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("ENCODE SOAP ", class(obj), "\n")
			}
			stop("Method not implemented yet!")
		}
)


################################################################################
# encoding functions

setMethod(f = "encodeXML",
		signature = signature(obj = "SosEventTime", sos = "SOS"),
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeXML]", class(obj), "\n")
			}
			
			.temporalOpsClass <- class(obj@temporalOps)
			if(!is.null(SosSupportedTemporalOperators()[[.temporalOpsClass]])) {
				.eventTime <- xmlNode(name = sosEventTimeName,
						namespace = sosNamespacePrefix)
				.temporalOpsXML <- encodeXML(obj = obj@temporalOps,
						sos = sos, verbose = verbose)
				.eventTime$children[[1]] <- .temporalOpsXML
				
				return(.eventTime)
			}
			else {
				stop(paste("temporalOps type not supported:",
								.temporalOpsClass))
			}
		}
)
setMethod(f = "encodeXML",
		signature = signature(obj = "SosEventTimeLatest", sos = "SOS"),
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeXML]", class(obj), "\n")
			}
			
			.eventTime <- xmlNode(name = sosEventTimeName,
					namespace = sosNamespacePrefix)
			.tmEquals <- xmlNode(name = ogcTempOpTMEqualsName,
					namespace = ogcNamespacePrefix)
			.propertyName <- xmlNode(name = ogcPropertyNameName,
					namespace = ogcNamespacePrefix)
			xmlValue(.propertyName) <- sosDefaultTempOpPropertyName
			.latestTime <- xmlNode(name = gmlTimeInstantName,
					namespace = gmlNamespacePrefix)
			.tpos <- xmlNode(name = gmlTimePositionName,
					namespace = gmlNamespacePrefix)
			xmlValue(.tpos) <- sosEventTimeLatestValue
			
			.latestTime$children[[1]] <- .tpos
			.tmEquals$children[[1]] <- .propertyName
			.tmEquals$children[[2]] <- .latestTime
			.eventTime$children[[1]] <- .tmEquals
			
			return(.eventTime)
		}
)

setMethod(f = "encodeXML",
		signature = signature(obj = "SosFeatureOfInterest", sos = "SOS"),
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeXML]", class(obj), "\n")
			}
			
			.foi <- xmlNode(name = sosFeatureOfInterestName,
					namespace = sosNamespacePrefix)
			
			# switch between objectIDs and spatialOps
			if(!any(is.na(obj@objectIDs))) {
				.ids <- lapply(X = obj@objectIDs, FUN = xmlNode,
						name = sosObjectIDName, namespace = sosNamespacePrefix)
				.foi <- addChildren(node = .foi, kids = .ids)
			}
			else if (!is.null(obj@spatialOps)) {
				.spOp <- encodeXML(obj = obj@spatialOps, sos = sos,
						verbose = verbose)
				.foi <- addChildren(node = .foi, kids = list(.spOp))
			}
			
			return(.foi)
		}
)

#
# to make just the time encoding interchangeable by users
#
setMethod(f = "encodeXML",
		signature = signature(obj = "POSIXt", sos = "SOS"),
		def = function(obj, sos, verbose = FALSE) {
			if(verbose) cat("[encodeXML] POSIXt with value", toString(obj),
						"\n")
			
			.formatted <- strftime(x = obj, format = sosTimeFormat(sos))
			
			if(verbose)
				cat("Formatted ", obj, " to ", .formatted)
			
			return(.formatted)
		}
)

#
# 
#
setMethod(f = "encodeKVP",
		signature = signature(obj = "SosEventTime", sos = "SOS"),
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("ENCODE KVP ", class(obj), "\n")
			}
			
			.temporalOpsKVP <- encodeKVP(obj = obj@temporalOps, sos = sos,
					verbose = verbose)
			return(.temporalOpsKVP)
		}
)

#
# 
#
setMethod(f = "encodeKVP",
		signature = signature(obj = "SosEventTimeLatest", sos = "SOS"),
		function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("ENCODE KVP ", class(obj), "\n")
			}
			# if eventTime is not given in GET binding, then the latest observation is returned
			return(NA_character_)
		}
)

#
# to make just the time encoding interchangeable by users
#
setMethod(f = "encodeKVP",
		signature = signature(obj = "POSIXt", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeKVP] POSIXt with value", toString(obj),
						"\n")
			
			.formatted <- strftime(x = obj, format = sosTimeFormat(sos))
			
			if(verbose)
				cat("Formatted ", obj, " to ", .formatted)
			
			return(.formatted)
		}
)


################################################################################
#
setMethod(f = "checkRequest",
		signature = signature(service = "SOS", operation = "SosDescribeSensor",
				verbose = "logical"),
		def = function(service, operation, verbose) {
			if(verbose) {
				cat("[checkRequest] Checking DescribeSensor... ")
			}
			
			# check if operation is for SOS and operation is DescribeSensor
			if(!(operation@service == sosService && 
						operation@request == sosDescribeSensorName)) {
				stop("Wrong input! Require classes 'SOS' as service and ''SosDescribeSensor' as operation.")
				return(FALSE)
			}
			
			# check if sensor in in listed in procedures
			.procedures = unique(unlist(sosProcedures(service)))
			.dsOperation <- sosOperation(service, sosDescribeSensorName)
			
			.procContained <- FALSE
			for (x in .procedures) {
				if(x == operation@procedure)
					.procContained <- TRUE
			}
			if(!.procContained)
				warning("Requested procedure ist not listed in capablities, service might return error!")
			
			
			# check if output format is supported by sos
			.oFSupported <- FALSE
			.supportedFormats <- .dsOperation@parameters[["outputFormat"]];
			.format <- gsub(operation@outputFormat, pattern = "\\&quot;",
					replacement = '"')
			
			if(!any(sapply(.supportedFormats,
							"==",
							.format),
					na.rm = TRUE)) {
				warning(paste("Outputformat has to be one of",
								paste(.supportedFormats, sep = ", ",
										collapse = " ")))
			}
			else {
				.oFSupported <- TRUE
			}
			
			# check if method is supported
			.methodSupported <- any(sapply(SosSupportedConnectionMethods(),
							"==", service@method))
			if(!.methodSupported)
				warning("Requested method type ist not listed in capablities for this operation, service might return error!")
			
			if(verbose) {
				cat("[checkRequest] Checks: procedure contained =",
						.procContained,
						", output supported =", .oFSupported,
						", method supported =", .methodSupported, "\n")
			}
			
			return(.procContained && .oFSupported && .methodSupported)
		})
		
setMethod(f = "checkRequest",
		signature = signature(service = "SOS", operation = "SosGetObservation",
				verbose = "logical"),
		def = function(service, operation, verbose) {
			# check if operation is for SOS and operation is DescribeSensor
			if(!(operation@service == sosService && 
						operation@request == sosGetObservationName)) {
				stop("Wrong input! Require classes 'SOS' as service and 'GetObservation' as operation.")
				return(FALSE)
			}
			
			# TODO implement checkRequest for GetObservation
			
			# check if given responseFormat is supported by the service
			
			# check if temporal operator and operand are a valid combination according to filter capabilities
			
			return(TRUE)
		}
)

setMethod(f = "checkRequest",
		signature = signature(service = "SOS",
				operation = "SosGetObservationById", verbose = "logical"),
		def = function(service, operation, verbose) {
			# check if operation is for SOS and operation is DescribeSensor
			if(!(operation@service == sosService && 
						operation@request == sosGetObservationByIdName)) {
				stop("Wrong input! Require classes 'SOS' as service and 'GetObservationById' as operation.")
				return(FALSE)
			}
			
			# TODO implement checkRequest for GetObservationById
			# see above!
			
			return(TRUE)
		}
)
