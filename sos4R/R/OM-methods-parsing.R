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
# Dispatch function for all exchangeable parsers for OM elements.
#
parseOM <- function(obj, sos, verbose = FALSE) {
	.om <- NULL
	
	# check if this is the outermost call and a document is given, not a node
	if(inherits(obj, xmlInternalDocumentName))
		.root <- xmlRoot(obj)
	else .root <- obj
	
	# switch submethods based on name
	.rootName <- xmlName(.root)
	
	.parsingFunction <- sosParsers(sos)[[.rootName]]
	if(!is.null(.parsingFunction)) {
		if(verbose) cat("[parseOM] rootName is", .rootName, "\n") #, "with: "); print(.parsingFunction)
		.om <- .parsingFunction(obj = .root, sos = sos, verbose = verbose)
		if(verbose) cat("[parseOM] Done!", .rootName, ":",
					#substr(toString(.om), 0, 200), "...\n")
					toString(.om), "\n")
	}
	else {
		warning(paste("[parseOM] No parsing function for given element", .rootName))
	}
	
	return(.om)
}

#
# Function extracts om:Obervation or om:Measurement from om:member.
#
parseObservationProperty <- function(obj, sos, verbose = FALSE) {
	# a member can only have one child element, so omit text node artefacts
	if(xmlSize(obj) >= 1) {
		.noneTexts <- .filterXmlChildren(obj, xmlTextNodeName, includeNamed = FALSE)
		.child <- .noneTexts[[1]]
		#.child <- xmlChildren(obj)[[1]]
		if(verbose) {
			cat("[parseObservationProperty] Parsing child of member:",
					xmlName(.child), "\n")
		}
		.mResult <- parseOM(.child, sos, verbose)
	}
	else {
		# no child, try href attribute
		if(verbose) cat("[parseObservationProperty] Member has no direct child!\n")
		
		.href <- xmlGetAttr(node = obj, name = "href", default = NA_character_)
		if(!is.na(.href)) {
			warning(paste("[parseObservationProperty] Only reference to Observation was returned:",
					.href))
			.mResult <- OmObservationProperty(href = .href)
		}
		else {
			warning("[parseObservationProperty] No Observation found in response!")
			.mResult <- OmObservationProperty()
		}
	}
	
	return(.mResult)
}

#
# om:Measurement
#
parseMeasurement <- function(obj, sos, verbose = FALSE) {
	if(verbose) cat("[parseMeasurement]\n")
	
	.samplingTime <- parseSamplingTime(obj = obj[[omSamplingTimeName]],
			format = sosTimeFormat(sos), verbose = verbose)
	
	# 52N SOS only returns om:Measurements (!) with procedure ids and observed 
	# properties in xlink:href
	.procedure <- xmlGetAttr(node = obj[[omProcedureName]], name = "href")
	.observedProperty <- SwePhenomenonProperty(
			href = xmlGetAttr(node = obj[[omObservedPropertyName]],
					name = "href"))
	
	.featureOfInterest <- parseFOI(obj[[omFeatureOfInterestName]], sos = sos,
			verbose = verbose)
	
	# must be GmlMeasure
	.result <- parseMeasure(obj[[omResultName]])
	
	# TODO optionals elements for OmMeasurement
	#.metadata
	#.resultTime
	#.resultQuality
	#.parameter
	
	.measurement <- OmMeasurement(samplingTime = .samplingTime,
			procedure = .procedure, observedProperty = .observedProperty,
			featureOfInterest = .featureOfInterest, result = .result)
	
	return(.measurement)
}

#
# om:Observation
#
parseObservation <- function(obj, sos, verbose = FALSE) {
	.id <- xmlGetAttr(node = obj, name = "id",
			default = NA_character_)
	if(verbose) cat("[parseObservation]", .id, "\n")
	
	# 52N SOS only returns om:Observation with procedure ids xlink:href
	.procedure <- xmlGetAttr(node = obj[[omProcedureName]], name = "href",
			default = NA_character_)
	
	.observedProperty <- parsePhenomenonProperty(obj[[omObservedPropertyName]],
			sos = sos, verbose = verbose)
	
	if(!is.null(obj[[omSamplingTimeName]])) {
		.samplingTime <- parseSamplingTime(obj = obj[[omSamplingTimeName]],
				format = sosTimeFormat(sos = sos), verbose = verbose)
	} else {
		warning("om:samplingTime is mandatory in om:Observation, but is missing!")
		.samplingTime <- NULL
	}
	
	if(!is.null(obj[[omFeatureOfInterestName]])) {
		.featureOfInterest <- parseFOI(obj[[omFeatureOfInterestName]],
				sos = sos, verbose = verbose)
	} else {
		warning("om:featureOfInterest is mandatory in om:Observation, but is missing!")
		.featureOfInterest <- NULL
	}
	
	# result parser is exchangeable
	.resultParsingFunction <- sosParsers(sos)[[omResultName]]
	.result <- .resultParsingFunction(obj[[omResultName]], sos, verbose)
	
	# optional elements
	if(!is.null(obj[[omResultTimeName]])) {
		.resultTime <- parseSamplingTime(obj = obj[[omResultTimeName]],
				format = sosTimeFormat(sos = sos), verbose = verbose)
	}
	else {
		.resultTime <- NULL
	}
	
	# TODO optionals elements for OmObservation
	#.metadata
	#.resultQuality
	#.parameter
	#.metadata
	
	.obs <- OmObservation(samplingTime = .samplingTime,
			procedure = .procedure, observedProperty = .observedProperty,
			featureOfInterest = .featureOfInterest, result = .result)
	
	return(.obs)
}

#
#
#
parseObservationCollection <- function(obj, sos, verbose) {
	# remove nodes other than member
	.members <- .filterXmlChildren(obj, omMemberName, includeNamed = TRUE)
	
	if(verbose) cat("[parseObservationCollection] with ", length(.members), 
				"element(s).\n")
	
	.env <- obj[[gmlBoundedByName]][[gmlEnvelopeName]]
	if(!is.null(.env)) {
		.boundedBy <- list(
				srsName = xmlGetAttr(.env, "srsName"),
				lowerCorner = xmlValue(.env[[gmlLowerCornerName]]),
				upperCorner = xmlValue(.env[[gmlUpperCornerName]]))
		
		if(verbose) cat("[parseObservationCollection] Parsed envelope:",
					toString(.boundedBy))
		
		if(sosSwitchCoordinates(sos)) {
			warning("Switching coordinates in envelope of ObservationCollection!")
			.origLC <- strsplit(x = .boundedBy[["lowerCorner"]], split = " ")
			.lC <- paste(.origLC[[1]][[2]], .origLC[[1]][[1]])
			.origUC <- strsplit(x = .boundedBy[["upperCorner"]], split = " ")
			.uC <- paste(.origUC[[1]][[2]], .origUC[[1]][[1]])
			.boundedBy <- list(srsName = xmlGetAttr(.env, "srsName"),
					lowerCorner = .lC, upperCorner = .uC)
		}
	}
	else {
		if(verbose) cat("[parseObservationCollection] Empty envelope! ")
		.boundedBy <- list()
	}
	
	.resultList <- lapply(.members, parseOM, sos, verbose)
	
	names(.resultList) <- lapply(.resultList, class)
	
	if(is.list(.resultList)) {
		.obsColl <- OmObservationCollection(members = .resultList,
				boundedBy = .boundedBy)
	}
	else {
		.obsColl <- OmObservationCollection(members = list(.resultList),
				boundedBy = .boundedBy)
	}
	
	if(verbose)
		cat("[parseObservationCollection] Done. Processed", length(.obsColl),
				"elements:", names(sosResult(.obsColl)), "\n")
	
	return(.obsColl)
}

#
# om:result
#
parseResult <- function(obj, sos, verbose = FALSE) {
	if(verbose) {
		cat("[parseResult]\n")
#		print(obj)
	}
	.result <- NULL
	
	.noneText <- .filterXmlChildren(node = obj, xmlTextNodeName,
			includeNamed = FALSE, verbose = verbose)
	
	if(verbose) {
		cat("[parseResult]", length(.noneText), " non-text nodes, names:",
				names(.noneText), "\n")
	}
	
	# Check if remaining element is there
	if(length(.noneText) == 0) {
		.children <- xmlChildren(obj)
		cat("[parseResult] No non-text nodes in result, returning NULL.\n")
		return(NULL)
	}
	
	# 52N SOS currently only returns swe:DataArrayDocument, but still check
	if(xmlName(.noneText[[1]]) != sweDataArrayName) {
		warning(paste("[parseResult] Parsing of given result is NOT supported:",
						xmlName(.noneText[[1]]), "-- only", sweDataArrayName,
						"can be parsed."))
	}
	else {
		if(verbose) cat("[parseResult] Parsing result with swe:DataArray.\n")
		
		# data array parser is exchangeable
		.dataArrayParsingFunction <- sosParsers(sos)[[sweDataArrayName]]
		.dataArray <- .noneText[[1]]
		.result <- .dataArrayParsingFunction(.dataArray, sos, verbose)
	}
	
	if(is.null(.result)) {
		stop("[parseResult] result is null! Given result:\n")
		print(obj)
	}
	else return(.result)
}


################################################################################
# not yet supported specializations (constraints):

parseGeometryObservation <- function(obj, sos, verbose = FALSE) {
	warning("Parsing of om:GeometryObservation is not implemented!")
	return(NA)
}

parseCategoryObservation <- function(obj, sos, verbose = FALSE) {
	warning("Parsing of om:CategoryObservation is not implemented!")
	return(NA)
}

parseCountObservation <- function(obj, sos, verbose = FALSE) {
	warning("Parsing of om:CountObservation is not implemented!")
	return(NA)
}

parseTruthObservation <- function(obj, sos, verbose = FALSE) {
	warning("Parsing of om:TruthObservation is not implemented!")
	return(NA)
}

parseTemporalObservation <- function(obj, sos, verbose = FALSE) {
	warning("Parsing of om:TemporalObservatio is not implemented!")
	return(NA)
}

parseComplexObservation <- function(obj, sos, verbose = FALSE) {
	warning("Parsing of om:ComplexObservation is not implemented!")
	return(NA)
}


################################################################################
# not exchangeable parsing functions:

#
# parse sos:featureOfInterest to according Element of GML or SA
#
parseFOI <- function(obj, sos, verbose = FALSE) {
	if(verbose) cat("[parseFOI] starting...\n")
	.foi <- NULL
	
	# has href attribute? if yes, use it!
	.href <- xmlGetAttr(node = obj, name = "href")
	if(!is.null(.href)) {
		if(verbose) cat("[parseFOI] referenced FOI:", .href, "\n")
		# feature is referenced
		.foi <- GmlFeatureProperty(href = .href)
	}
	else {		
		# feature is available in the element
		.noneTexts <- .filterXmlChildren(obj, xmlTextNodeName,
				includeNamed = FALSE)
		.feature <- .noneTexts[[1]]
		.name <- xmlName(.feature)
		
		if(verbose) cat("[parseFOI] inline FOI:", .name, "\n")
		
		if(.name == saSamplingPointName) {
			.sp <- parseSamplingPoint(.feature, sos = sos)
			.foi <- GmlFeatureProperty(feature = .sp)
		}
		else if (.name == saSamplingSurface) {
			warning("[parseFOI] No parsing for sa:SamplingSurface implemented!")
			.foi <- GmlFeatureProperty(href = .name)
		}
		else if (.name == gmlFeatureCollectionName) {
			.foi <- parseFeatureCollection(.feature, sos = sos)
		}
		else {
			warning("[parseFOI] No parsing for given feature implemented!")
			.foi <- GmlFeatureProperty(href = .name)
		}
	}
	
	return(.foi)
}

#
# create according GmlTimeObject from om:samplingTime
#
parseSamplingTime <- function(obj, format, verbose = FALSE) {
	if(verbose) cat("[parseSamplingTime]\n")
	
	.tiXML <- xmlChildren(obj)[[gmlTimeInstantName]]
	.tpXML <- xmlChildren(obj)[[gmlTimePeriodName]]
	.timeObject <- NULL
	if(!is.null(.tiXML)) {
		if(verbose) cat("[parseSamplingTime] time instant.\n")
		.timeObject <- parseTimeInstant(obj = .tiXML, format = format)
	}
	else if(!is.null(.tpXML)) {
		if(verbose) cat("[parseSamplingTime] time period.\n")
		.timeObject <- parseTimePeriod(obj = .tpXML, format = format)
	}
	else {
		warning(paste("Could not create GmlTimeObject from given samplingTime,", 
					" require gml:TimeInstant or gml:TimePeriod as children."))
		.timeObject <- GmlTimeInstant(timePosition = GmlTimePosition(
						time = as.POSIXct(x = NA)))
	}
	
	return(.timeObject)
}

