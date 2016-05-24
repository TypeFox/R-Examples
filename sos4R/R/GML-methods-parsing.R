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
# Created: 2010-09-15                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################


#
#
#
parsePosition <- function(obj, sos) {
	.position <- NULL
	
	# has href attribute?
	.href <- xmlGetAttr(node = obj, name = "href")
	if(!is.null(.href)) {
		# position is referenced
		.position <- GmlPointProperty(href = .href)
	}
	else {
		# must be point
		.position <- GmlPointProperty(point = parsePoint(obj[[gmlPointName]],
						sos = sos))
	}
	
	return(.position)
}

#
#
#
parsePoint <- function(obj, sos) {
	.point <- NA
	.pos <- obj[[gmlPosName]]
	
	.posString <- xmlValue(.pos)
	
	if(sosSwitchCoordinates(sos)) {
		warning("Switching coordinates in Point!")
		.orig <- strsplit(x = .posString, split = " ")
		.posString <- paste(.orig[[1]][[2]], .orig[[1]][[1]])
	}
	
	# optional attributes:
	.srsName <- xmlGetAttr(node = .pos, name = "srsName",
			default = NA_character_)
	.srsDimension <- xmlGetAttr(node = .pos, name = "srsDimension",
			default = NA_integer_)
	.axisLabels <- xmlGetAttr(node = .pos, name = "axisLabels",
			default = NA_character_)
	.uomLabels <- xmlGetAttr(node = .pos, name = "uomLabels",
			default = NA_character_)
	
	.pos <- GmlDirectPosition(pos = .posString, srsName = .srsName,
			srsDimension = .srsDimension, axisLabels = .axisLabels,
			uomLabels = .uomLabels)
	.point <- GmlPoint(pos = .pos)	
	
	return(.point)
}

#
#
#
parseTimeInstant <- function(obj, format) {
	.timePosXML <- .filterXmlChildren(node = obj, 
			childrenName = gmlTimePositionName, includeNamed = TRUE)[[1]]
	
	.timePos <- parseTimePosition(obj = .timePosXML,
			format = format)
	
	#optionals
	.id = xmlGetAttr(node = obj, name = "id",
			default = NA_character_)
	.frame = xmlGetAttr(node = obj, name = "frame",
			default = as.character(NA))
	.noneTexts <- .filterXmlChildren(node = obj, gmlRelatedTimeName)
	if(!is.null(.noneTexts))
		.relatedTimes <- .noneTexts
	else
		.relatedTimes = list()
	
	.ti <- GmlTimeInstant(timePosition = .timePos, id = .id,
			relatedTimes = .relatedTimes, frame = .frame)
	return(.ti)
}

#
#
#
parseTimeInstantProperty <- function(obj, format) {
	.timeProp <- NULL
	
	# check if reference or inline phenomenon
	.href <- xmlGetAttr(node = obj, name = "href")
	if(!is.null(.href)) {
		.timeProp <- GmlTimeInstantProperty(href = .href)
	}
	else {
		.noneText <- .filterXmlChildren(node = obj, xmlTextNodeName,
				includeNamed = FALSE)
		.time <- parseTimeInstant(obj = .noneText[[1]], format = format)
		.timeProp <- GmlTimeInstantProperty(time = .time)
	}
	
	return(.timeProp)
}

#
# 
#
parseTimePosition <- function(obj, format) {
	.time <- strptime(xmlValue(obj), format)
	
	# optional:
	.frame <- xmlGetAttr(node = obj, name = "frame", default = NA_character_)
	.calendarEraName <- xmlGetAttr(node = obj, name = "calendarEraName",
			default = NA_character_)
	.indeterminatePosition <- xmlGetAttr(node = obj,
			name = "indeterminatePosition", default = NA_character_)
	
	.timePosition <- GmlTimePosition(time = .time, frame = .frame,
			calendarEraName = .calendarEraName,
			indeterminatePosition = .indeterminatePosition)
}

#
#
#
parseTimePeriod <- function(obj, format) {
	.timeObject <- NULL
	
	# optionals
	.id = xmlGetAttr(node = obj, name = "id",
			default = NA_character_)
	.frame = xmlGetAttr(node = obj, name = "frame",
			default = as.character(NA))
	.noneTexts <- .filterXmlChildren(node = obj, gmlRelatedTimeName)
	if(!is.null(.noneTexts))
		.relatedTimes <- .noneTexts
	else
		.relatedTimes = list()
	
	# TODO parse gml:timeLength
	.duration <- NA_character_
	.timeInterval <- NULL
	
	# begin and end
	if(!is.null(obj[[gmlBeginName]]) || !is.null(obj[[gmlEndName]])) {
		.begin <- parseTimeInstantProperty(obj = obj[[gmlBeginName]],
				format = format)
		.end <- parseTimeInstantProperty(obj[[gmlEndName]], format = format)
		
		.timeObject <- GmlTimePeriod(begin = .begin, end = .end, duration = .duration,
				timeInterval = .timeInterval, id = .id,
				relatedTimes = .relatedTimes, frame = .frame)
	}
	# beginPosition and endPosition
	else if(!is.null(obj[[gmlBeginPositionName]])
			|| !is.null(obj[[gmlEndPositionName]])) {
		.beginPosition <- parseTimePosition(
				obj = obj[[gmlBeginPositionName]],
				format = format)
		.endPosition <- parseTimePosition(
				obj = obj[[gmlEndPositionName]],
				format = format)
		
		.timeObject <- GmlTimePeriod(beginPosition = .beginPosition,
				endPosition = .endPosition, duration = .duration,
				timeInterval = .timeInterval, id = .id,
				relatedTimes = .relatedTimes, frame = .frame)
	}
	
	return(.timeObject)
}

#
#
#
parseTimeGeometricPrimitiveFromParent <- function(obj, format) {
	.tiXML <- xmlChildren(obj)[[gmlTimeInstantName]]
	.tpXML <- xmlChildren(obj)[[gmlTimePeriodName]]
	.timeObject <- NULL
	if(!is.null(.tiXML)) {
		.timeObject <- parseTimeInstant(obj = .tiXML, format = format)
	}
	else if(!is.null(.tpXML)) {
		.timeObject <- parseTimePeriod(obj = .tpXML, format = format)
	}
	else {
#		warning(paste("Could not create time from given samplingTime,", 
#						" require gml:TimeInstant or gml:TimePeriod as children."))
		.timeObject <- GmlTimeInstant(timePosition = GmlTimePosition(
						time = as.POSIXct(x = NA)))
	}
	
	return(.timeObject)
}

#
#
#
parseFeatureCollection <- function(obj, sos) {
	.members <- .filterXmlChildren(node = obj,
			childrenName = gmlFeatureMemberName)
	
	.id <- xmlGetAttr(node = obj, name = "id", default = NA_character_)
	
	if(length(.members) > 0) {
		.members <- lapply(.members, .parseFeatureMember, sos = sos)
		
		.fc <- GmlFeatureCollection(featureMembers = .members, id = .id)
		return(.fc)
	}
	else {
		warning("gml:FeatureCollection could not be parsed, only gml:featureMember elements can be handled.")
		return(NULL)
	}
}
.parseFeatureMember <- function(obj, sos) {
	.noneTexts <- .filterXmlOnlyNoneTexts(obj)
	.member <- .noneTexts[[1]]
	
	.name <- xmlName(.member)
	
	if(.name == saSamplingPointName) {
		.sp <- parseSamplingPoint(.member, sos = sos)
		.member.parsed <- GmlFeatureProperty(feature = .sp)
	}
	else if (.name == gmlFeatureCollectionName) {
		.member.parsed <- parseFeatureCollection(.member, sos = sos)
	}
	else {
		warning("No handling for given gml:featureMember available, only sa:SamplingPoint is supported!")
	}
	return(.member.parsed)
}

#
#
#
parseMeasure <- function(obj) {
	.value <- as.numeric(xmlValue(obj))
	.uom <- xmlGetAttr(node = obj, name = "uom", default = NA_character_)
	.result <- GmlMeasure(.value, .uom)
	return(.result)
}

