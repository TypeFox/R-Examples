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
# Created: 2010-09-21                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

################################################################################
# read function, a convenience function using non-SWE names for stuff
#
# parameters:
#		sos:
#				the sos to query
# 		time = c():	
#				either one time as POSIXt (or character which is tried to be
#				parsed to POSIXt) and forms a time instant, or two times
#				(character or POSIXt) that form an interval
#		mergeResult = FALSE:
#				boolean flag that turns on merging of the result into one
#				dataframe, the single observations cannot be accesses,
#				essentially only gives all [measurement,observation]@result
#		addLocation = FALSE
#				trigger if the location from the  should be added to the result
#				table fom the featureOfInterest
#
# * Function does not include "offering", which is/are instead automatically
#	detected... but what if one sensor is in several offerings?
#
# * Function returns data.frame, NOT OmObservation or OmMeasurement classes
#
# * if querying several procedures with different positions, put in them into 
#	one data.frame, or a list?
#
read.sos <- function(sos,
		sensors = NA_character_,
		phenomena = NA_character_,
		bbox = NA_character_, # one, or several?
		times = NA_character_, # one, or several?
		mergeResult = FALSE,
		addLocation = FALSE,
		verbose = FALSE) {
	warning("Method is not implemented yet!")
}


# TODO method creates a list (whose items can are named and can be used to read data,
# or just used in function read.sos)
# for all available and valid combinations of off/foi/phen/proc
#
#sosTimeSeries <- function(sos) {
#	warning("Method is not implemented yet!")
#}


# TODO get all the matching ungiven parameters for that are available for a set of 
# parameters that are given, e.g. all offerings that offer observed property
# "A" for feature of interest "X", or all procedures measuring an observed
# property "B".
#
#sosMatching <- function(offering = NA, procedure = NA, observedProperty = NA,
#		foi = NA) {
#	warning("Method is not implemented yet!")
#}


################################################################################
# conversion methods
#
sosConvertTime <- function(x, sos) {
	.t <- as.POSIXct(x = strptime(x = x, format = sosTimeFormat(sos = sos)))
	return(.t)
}
sosConvertDouble <- function(x, sos) {
	return(as.double(x = x))
}
sosConvertString <- function(x, sos) {
	return(as.character(x = x))
}
sosConvertLogical <- function(x, sos) {
	return(as.logical(x = x))
}


################################################################################
# convenience functions
#
setMethod(f = "sosCreateTimeInstant",
		signature = signature(sos = "SOS", time = "POSIXt"),
		def = function(sos, time, frame, calendarEraName,
				indeterminatePosition) {
#			.time <- format(time, sosTimeFormat(sos))
			.timePos <- GmlTimePosition(
#					time = strptime(.time, sosTimeFormat(sos)),
					time = time,
					frame = frame, calendarEraName = calendarEraName,
					indeterminatePosition = indeterminatePosition)
			.ti <- GmlTimeInstant(timePosition = .timePos)
			return(.ti)
		}
)

#
#
#
setMethod(f = "sosCreateTimePeriod",
		signature = signature(sos = "SOS", begin = "POSIXt", end = "POSIXt"),
		def = function(sos, begin, end, frame, calendarEraName,
				indeterminatePosition, duration, timeInterval) {
#			.tf <- sosTimeFormat(sos)
			.beginPos <- GmlTimePosition(
#					time = strptime(format(begin, .tf), .tf),
					time = begin,
					frame = frame, calendarEraName = calendarEraName,
					indeterminatePosition = indeterminatePosition
			)
			.endPos <- GmlTimePosition(
#					time = strptime(format(end, .tf), .tf),
					time = end,
					frame = frame, calendarEraName = calendarEraName,
					indeterminatePosition = indeterminatePosition
			)
			.tp <- GmlTimePeriod(beginPosition = .beginPos,
					endPosition = .endPos, duration = duration,
					timeInterval = timeInterval)
			return(.tp)
		}
)

#
#
#
setMethod(f = "sosCreateEventTimeList",
		signature = signature(time = "GmlTimeGeometricPrimitive"),
		def = function(time, operator) {
			.et <- list(sosCreateEventTime(time = time, operator = operator))
			return(.et)
		}
)

#
#
#
setMethod(f = "sosCreateTime",
		signature = signature(sos = "SOS", time = "character"),
		def = function(sos, time, operator) {
			if(regexpr(pattern = "::", text = time) > -1) {
				.l <- .sosCreateEventTimeListFromPeriod(sos = sos, time = time,
						operator = operator, seperator = "::")
			}
			else if(regexpr(pattern = "P", text = time) > -1) {
				.l <- .sosCreateEventTimeListFromISOPeriod(sos = sos, 
						time = time, operator = operator)
			}
			else if(regexpr(pattern = "/", text = time) > -1) {
				.l <- .sosCreateEventTimeListFromPeriod(sos = sos, time = time,
						operator = operator, seperator = "/")
			}

			return(.l)
		}
)

.sosCreateEventTimeListFromPeriod <- function(sos, time, operator, seperator) {
	.times <- strsplit(x = time, split = seperator)[[1]]
	.start <- .times[[1]]
	if(length(.times) > 1)
		.end <- .times[[2]]
	else .end <- NULL
	
#	print(.start); print(.end);	print(nchar(.start)); print(nchar(.end));
#	str(.start); print(.end);
	
	if(is.null(.start) && is.null(.end)) {
		warning("Both start and endtime are null based on given time. Returning empty list!")
		return(list())
	}
	else {
		if(is.null(.end)) {
			# no end time:
			.ti <- sosCreateTimeInstant(sos = sos, time = as.POSIXct(.start))
			.l <- sosCreateEventTimeList(time = .ti,
					operator = SosSupportedTemporalOperators()[[ogcTempOpTMAfterName]])
		}
		else if(nchar(.start) > 0) {
			.tp <- sosCreateTimePeriod(sos = sos, begin = as.POSIXct(.start),
					end = as.POSIXct(.end))
			.l <- sosCreateEventTimeList(.tp)
		}
		else if(nchar(.start) < 1) {
			# no start time:
			.ti <- sosCreateTimeInstant(sos = sos, time = as.POSIXct(.end))
			.l <- sosCreateEventTimeList(time = .ti,
					operator = SosSupportedTemporalOperators()[[ogcTempOpTMBeforeName]])
		}
	}
	
	return(.l)
}

.sosCreateEventTimeListFromISOPeriod <- function(sos, time, operator) {
#	* 2005-08-09T18:31:42P3Y6M4DT12H30M17S: bestimmt eine Zeitspanne von 3 Jahren, 6 Monaten, 4 Tagen 12 Stunden, 30 Minuten und 17 Sekunden ab dem 9. August 2005 "kurz nach halb sieben Abends".
#	* P1D: "Bis morgen zur jetzigen Uhrzeit." Es koennte auch "PT24H" verwendet werden, doch erstens waeren es zwei Zeichen mehr, und zweitens wuerde es bei der Zeitumstellung nicht mehr zutreffen.
#	* P0003-06-04T12:30:17
#	* P3Y6M4DT12H30M17S: gleichbedeutend mit dem ersten Beispiel, allerdings ohne ein bestimmtes Startdatum zu definieren
#	* PT72H: "Bis in 72 Stunden ab jetzt."
#	* 2005-08-09P14W: "Die 14 Wochen nach dem 9. August 2005."
#	* 2005-08-09/2005-08-30
#	* 2005-08-09--30
#	* 2005-08-09/30: "Vom 9. bis 30. August 2005."
	
	warning("Function .sosCreateEventTimeListFromISOPeriod not implemented yet!")
}

#
#
#
setMethod(f = "sosCreateEventTime",
		signature = signature(time = "GmlTimeGeometricPrimitive"),
		def = function(time, operator) {
			
			if(operator == ogcTempOpTMAfterName) {
				.tOps <- TM_After(time = time)
			}
			else if(operator == ogcTempOpTMBeforeName) {
				.tOps <- TM_Before(time = time)
			}
			else if(operator == ogcTempOpTMDuringName) {
				.tOps <- TM_During(time = time)
			}
			else if(operator == ogcTempOpTMEqualsName) {
				.tOps <- TM_Equals(time = time)
			}
			else {
				stop(paste("Given operator", operator, "is not supported,",
								"choose one of",
								toString(SosSupportedTemporalOperators())))
			}
			
			.et <- SosEventTime(.tOps)
			return(.et)
		}
)

#
#
#
setMethod(f = "sosCreateFeatureOfInterest",
		signature = signature(),
		def = function(objectIDs, spatialOps, bbox, srsName) {
			# switch cases, either objectIDs or one of the spatialOps shortcuts
			if(!any(is.na(objectIDs))) {
				.foi <- SosFeatureOfInterest(objectIDs = objectIDs)
			}
			else if (!is.null(spatialOps)) {
				.foi <- SosFeatureOfInterest(spatialOps = spatialOps)
			}
			else if(!is.null(bbox)) {
				if(is.matrix(bbox)) {
					.env <- GmlEnvelope(
							lowerCorner = GmlDirectPositionLatLon(lat = bbox[2,1],
									lon = bbox[1,1]),
							upperCorner = GmlDirectPositionLatLon(lat = bbox[2,2],
									lon = bbox[1,2]),
							srsName = srsName)
					.bbox <- OgcBBOX(envelope = .env)
					.foi <- SosFeatureOfInterest(spatialOps = .bbox)
				}
				else {
					stop("bbox must be matrix!")
				}
			}
			else {
				stop("At least one of objectIDs or spatialOps has to be set!")
			}
			
			return(.foi)
		}
)
		
#
#
#
setMethod(f = "sosCreateBBOX",
		signature = signature(lowLat = "numeric", lowLon = "numeric",
				uppLat = "numeric", uppLon = "numeric"),
		def = function(lowLat, lowLon, uppLat, uppLon, srsName,
			srsDimension = NA_integer_, axisLabels = NA_character_,
			uomLabels = NA_character_,
			propertyName = sosDefaultSpatialOpPropertyName) {
		.env <- GmlEnvelope(
				lowerCorner = GmlDirectPosition(
						pos = paste(lowLat, lowLon, sep = " ")),
				upperCorner = GmlDirectPosition(
						pos = paste(uppLat, uppLon, sep = " ")),
				srsName = srsName, srsDimension = srsDimension,
				axisLabels = axisLabels, uomLabels = uomLabels)
		
		.bbox <- OgcBBOX(propertyName = propertyName, envelope = .env)
		return(.bbox)
		}
)

#
#
#
setMethod(f = "sosCreateBBoxMatrix",
		signature = signature(lowLat = "numeric", lowLon = "numeric",
				uppLat = "numeric", uppLon = "numeric"),
		def = function(lowLat, lowLon, uppLat, uppLon) {
			.m <- matrix(data = c(lowLon, lowLat, uppLon, uppLat),
					nrow = 2, ncol = 2,
					dimnames = list(
							c("longitude", "latitude"),
							c("lowerCorner", "upperCorner")))
			return(.m)
		}
)

################################################################################
# MISC

#
#
#
setMethod(f = "sosCapabilitiesDocumentOriginal",
		signature = signature(sos = "SOS"),
		def = function(sos) {
			.gc <- OwsGetCapabilities(service = sosService,
					acceptVersions = c(sos@version))
			.responseString = sosRequest(sos = sos, request = .gc,
					verbose = sos@verboseOutput, inspect = FALSE)
			.response <- xmlParseDoc(.responseString, asText = TRUE)
			return(.response)
		}
)

#
# Helper function to get the capabilities URL, e.g. in Sweave documents
#
setMethod(f = "sosCapabilitiesUrl",
		signature = signature(sos = "SOS"),
		def = function(sos) {
			.gc <- OwsGetCapabilities(service = sosService,
					acceptVersions = c(sos@version))
			
			.request <- paste0(sosUrl(sos), "?", encodeRequestKVP(.gc, sos))
			return(.request)
		}
)

#
# helper methods for exception response handling
#
.isExceptionReport <- function(document) {
	if(owsExceptionReportName == xmlName(xmlRoot(document)))
		return(TRUE)
	else
		return(FALSE)
}
.handleExceptionReport <- function(sos, obj) {
	if(sos@verboseOutput) warning("Received ExceptionReport!")
	
	.parsingFunction <- sosParsers(sos)[[owsExceptionReportName]]
	.er <- .parsingFunction(obj)
	if(class(.er) == "OwsExceptionReport")
		warning(toString(.er))
	return(.er)
}

#
#
#
.createLatestEventTime <- function(verbose = FALSE) {
	if(verbose) cat("Creating non-standard event time 'latest'\n")
	.et <- SosEventTimeLatest()
	return(.et)
}

#
# encoding functions that just pass given content along...
#
setMethod(f = "encodeXML", signature = signature(obj = "XMLNode", sos = "SOS"),
		def = function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeXML] from XMLNode\n")
			}
			return(obj)
		}
)
setMethod(f = "encodeXML", signature = signature(obj = "XMLInternalElementNode", sos = "SOS"),
		def = function(obj, sos, verbose = FALSE) {
			if(verbose) {
				cat("[encodeXML] from XMLInternalElementNode: just returning it.\n")
			}
			return(obj)
		}
)
setMethod(f = "encodeXML", signature = signature(obj = "character", sos = "SOS"),
		def = function(obj, sos, addNamespaces = FALSE, verbose = FALSE) {
			if(verbose) cat("[encodeXML] from character string\n")
			
			if(isXMLString(obj)) {
				#FIXME this just won't work, see testing.R, section "encode xml character string (again)"
				if(addNamespaces) {
					if(verbose) cat("[encodeXML] Namespace hack for character string, trying to replace 'result>'!\n")
					.hack <- 'result xmlns:sos="http://www.opengis.net/sos/1.0" xmlns:xsi="http://www.w3.org/2001/XMLSchema-instance" xmlns:ows="http://www.opengis.net/ows/1.1" xmlns:om="http://www.opengis.net/om/1.0" xmlns:ogc="http://www.opengis.net/ogc" xmlns:gml="http://www.opengis.net/gml">'
					.hackedString <- sub(pattern = "result>",
							replacement = .hack,
							x = obj)
					.xml <- xmlTreeParse(.hackedString, asText = TRUE,
							useInternalNodes = FALSE)
				}
				else {
					.xml <- xmlParseString(obj)
				}
								
				if(verbose) {
					cat("[encodeXML] Created XML from string:\n", toString(.xml))
				}
				return(.xml)
			}
			else {
				warning(paste("[encodeXML] Could not encode given character string as XML!",
								" Character string: '", obj, "'", sep = ""))
			}
		}
)

################################################################################
# Helper functions for OWS exceptions, e.g. to get the meaning of an exception
# code.
#
setMethod(f = "sosExceptionCodeMeaning",
		signature = c(exceptionCode = "character"),
		def = function(exceptionCode) {
			.meaning <- as.character(
					.owsStandardExceptions[
							.owsStandardExceptions$exceptionCode==exceptionCode,
							2])
			return(.meaning)
		}
)

################################################################################
# get coordiante refernce system CRS
#
#
setMethod(f = "sosGetCRS",
		signature = c(obj = "character"),
		def = function(obj, verbose = FALSE) {
			if(verbose) cat("[sosGetCRS] from '", obj, "'\n", sep = "")
			
			# get the position of EPSG code
			.split <- strsplit(as.character(obj), split = ":")
			.idx <- which(toupper(.split[[1]]) == "EPSG")
			if(length(.idx) == 0) {
				# possibly versioned, try one index higher?
				warning(paste("Could not create CRS from the given object:", obj))
				return(NULL)
			}
			.epsg <- .split[[1]][[length(.split[[1]])]]
			
			.initString <- paste("+init=epsg", .epsg, sep = ":")
			
			if(verbose) cat("[sosGetCRS] .initString:", .initString, "\n")
			
			.rgdal <- require("rgdal", quietly = TRUE)
			if(!.rgdal)
				# if(!("rgdal" %in% .packages())) does only check loaded pkgs
				warning("[sosGetCRS] rgdal not present: CRS values will not be validated.",
						immediate. = TRUE)
			else
				if(verbose) cat("[sosGetCRS] rgdal loaded! \n")
				
			.crs <- NULL
			tryCatch({
						.crs <- CRS(.initString)
					}, error = function(err) {
						warning("[sosGetCRS] error was detected, probably the ",
								"EPSG code ", .epsg, " is not recognized ", 
								"(returning NULL):", toString(err))
					})
			
			if(verbose) {
				cat("[sosGetCRS] found: ")
				show(.crs)
			}
			
			return(.crs)
		}
)
setMethod(f = "sosGetCRS",
		signature = c(obj = "OmObservationCollection"),
		def = function(obj, verbose = FALSE) {
			.l <- lapply(X = obj, FUN = sosGetCRS, verbose = verbose)
			.l <- unique(.l)
			
			if(length(.l) == 1)
				return(.l[[1]])
			else return(.l)
		}
)
setMethod(f = "sosGetCRS",
		signature = c(obj = "OmObservation"),
		def = function(obj, verbose = FALSE) {
			.crs <- .getCRSfromOM(obj)
			return(.crs)
		}
)
setMethod(f = "sosGetCRS",
		signature = c(obj = "OmMeasurement"),
		def = function(obj, verbose = FALSE) {
			.crs <- .getCRSfromOM(obj)
			return(.crs)
		}
)
setMethod(f = "sosGetCRS",
		signature = c(obj = "SosObservationOffering"),
		def = function(obj, verbose = FALSE) {
			.srsName <- sosBoundedBy(obj)[["srsName"]]
			if(is.null(.srsName))
				.crs <- NULL
			else .crs <- sosGetCRS(.srsName, verbose = verbose)
			return(.crs)
		}
)
setMethod(f = "sosGetCRS",
		signature = c(obj = "SOS"),
		def = function(obj, verbose = FALSE) {
			.offs <- sosOfferings(obj)
			.crss <- lapply(.offs, sosGetCRS, verbose = verbose)
			if(length(.crss) == 1)
				return(.crss[[1]])
			return(.crss)
		}
)
setMethod(f = "sosGetCRS",
		signature = c(obj = "list"),
		def = function(obj, verbose = FALSE) {
			.crs <- lapply(X = obj, FUN = sosGetCRS, verbose = verbose)
			return(.crs)
		}
)

.getCRSfromOM <- function(obj) {
	.char <- as.vector(sosCoordinates(obj)[[sosDefaultColumnNameSRS]])
	.l <- sapply(X = .char, FUN = sosGetCRS)
	.l <- unique(.l)
	
	if(length(.l) == 1)
		return(.l[[1]])
	else return(.l)
}


################################################################################
#
# ", *, :, /, <, >, ?, \, and |
#
.cleanupFileName <- function(obj) {
	.clean <- gsub(
			pattern = "[\\/:\"|?<>*,]+",
			x = obj,
			replacement = "_")
	return(.clean)
}

.illegalColumnNameCharacters <- list("\\[", "\\]", "@", "\\$", "~",
		"\\+", "-", "\\*")
.illegalColumnNameEscapeCharacter <- "."

.cleanupColumnName <- function(name) {
	# replace illegal characters
	.name <- name
	
	for (.x in .illegalColumnNameCharacters) {
		# replace multiple escape characters with one
		.name <- gsub(pattern = .x,
				replacement = .illegalColumnNameEscapeCharacter,
				x = .name)
	}
	
	.name <- gsub(pattern = paste("(\\",
					.illegalColumnNameEscapeCharacter, ")+", sep = ""),
			replacement = .illegalColumnNameEscapeCharacter, x = .name)
	return(.name)
}


################################################################################
# Access to CHANGES and NEWS shipping with the package
#
#
sosChanges <- function() {
	.path <- paste(find.package("sos4R", lib.loc = NULL), "CHANGES",
			sep = "\\")
	.con <- file(.path)
	.lines <- readLines(.con)
	close(.con)	
	
	cat(.lines, sep = "\n")
}

#
#
#
sosNews <- function() {
	.path <- paste(find.package("sos4R", lib.loc = NULL), "NEWS",
			sep = "\\")
	.con <- file(.path)
	.lines <- readLines(.con)
	close(.con)	
	
	cat(.lines, sep = "\n")
}

#
# based on vignette-function
#
sosCheatSheet <- function() {
	.path <- paste(find.package("sos4R", lib.loc = NULL), "doc",
			.sosCheatSheetDocumentName, sep = "\\")
	
	.z <- list(file = .sosCheatSheetDocumentName, pdf = .path)
	.z$topic <- "sos4R Cheat Sheet"
	class(.z) <- "vignette"

	return(.z)
}
