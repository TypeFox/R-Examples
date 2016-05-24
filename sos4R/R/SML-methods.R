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
SensorML <- function(xml, coords = data.frame(), id = NA_character_,
		name = NA_character_, description = NA_character_,
		boundedBy = matrix()) {
	new("SensorML", xml = xml, coords = coords, id = id, name = name,
			description = description, boundedBy = boundedBy)
}

.xPathToken <- "@@@"
.smlXPathIdentifier <- paste(
		"//sml:System/sml:identification//sml:identifier/sml:Term[@definition='urn:ogc:def:identifier:OGC::",
		.xPathToken, 
		"' or @definition='urn:ogc:def:identifier:OGC:1.0:",
		.xPathToken,
		"' or @definition='urn:ogc:def:identifier:OGC:", # technically this is invalid, but common
		.xPathToken,
		"']/sml:value/text()", sep = "")
.smlXPathDescription <- "//sml:System/gml:description/text()"
.smlXPathPosition <- "//sml:System/sml:position/swe:Position"
.smlXPathObservedBBox <- "//swe:field[@name='observedBBOX']/swe:Envelope"


#
# parseSensorML(mySensor@xml, sos = mySOS, verbose = TRUE)
#
parseSensorML <- function(obj, sos, verbose = FALSE) {
	.root <- xmlRoot(obj)
	if(verbose) cat("[parseSensorML] Starting... \n")
	
	.id <- .smlIdentifier(.root, "uniqueID", verbose = verbose)
	.shortName <- .smlIdentifier(.root, "shortName", verbose = verbose)
	.descrNodeSet <- getNodeSet(doc = .root, path = .smlXPathDescription,
			namespaces = .sosNamespaceDefinitionsSML)
	if(is.null(.descrNodeSet))
		.description <- NA_character_
	else
		.description <- xmlValue(.descrNodeSet[[1]])
	if(verbose) cat("[parseSensorML] Got ID", .id, "and shortName", .shortName,
				"and description", .description, "\n")

	# bounded by
	if(verbose) cat("[parseSensorML] Parsing boundedBy from",
				.smlXPathObservedBBox, "\n")
	.observedBBox <- getNodeSet(doc =.root,
			path = .smlXPathObservedBBox,
			namespaces = .sosNamespaceDefinitionsSML)[[1]]
	if(!is.null(.observedBBox)) {
		.referenceFrame <- .observedBBox[["referenceFrame"]]
		.llVector <- parseVector(.observedBBox[["lowerCorner"]][["Vector"]],
				sos = sos, verbose = verbose)
		.uuVector <- parseVector(.observedBBox[["upperCorner"]][["Vector"]],
				sos = sos, verbose = verbose)
		.bb <- matrix(c(.llVector[["x"]][["value"]],
						.llVector[["y"]][["value"]],
						.uuVector[["x"]][["value"]],
						.uuVector[["y"]][["value"]]),
				ncol = 2,
				dimnames = list(c("coords.lon", "coords.lat"),
						c("min", "max")))
		.oldAttrs <- attributes(.bb)
		attributes(.bb) <- c(.oldAttrs,
				list(referenceFrame = .referenceFrame))
		
		if(verbose) cat("[parseSensorML] Parsed bounding box: ", toString(.bb),
							"\n")
	}
	else {
		.bb <- matrix()
		if(verbose) cat("[parseSensorML] No boundedBy element found, bbox is ",
					.bb, "\n")
	}
	
	# coordinates
	if(verbose) cat("[parseSensorML] Parsing coordinates from",
				.smlXPathPosition, "\n")
	.xmlPosition <- getNodeSet(doc = .root, path = .smlXPathPosition,
			namespaces = .sosNamespaceDefinitionsSML)
	if(length((.xmlPosition)) > 0) {
		.xmlPosition <- .xmlPosition[[1]]
		.position <- parseSwePosition(.xmlPosition, sos = sos,
				verbose = verbose)
		.referenceFrame = attributes(.position)[["referenceFrame"]]			
		.uom <- lapply(.position, "[[", "uomCode")
		names(.uom) <- lapply(.position, "[[", "axisID")
		.name <- lapply(.position, "[[", "name")
		names(.name) <- lapply(.position, "[[", "axisID")
		
		.values <- lapply(.position, "[[", "value")
		names(.values) <- lapply(.position, "[[", "axisID")
		if(any(is.na(names(.values)))) {
			warning("[parseSensorML] No axisID given, cannot name data.frame with them, trying 'name'.")
			names(.values) <- lapply(.position, "[[", "name")
		}
		
		if(verbose) {
			cat("[parseSensorML] names: ", names(.values), "\n")
			cat("[parseSensorML] values: ", toString(.values),	"\n")
		}
		
		.coords <- data.frame(.values)
		.oldAttrs <- attributes(.coords)
		attributes(.coords) <- c(as.list(.oldAttrs),
				list(referenceFrame = .referenceFrame,
						uom = .uom, name = .name))
		
		if(!is.na(.id))
			row.names(.coords) <- .id
		if(verbose) cat("[parseSensorML]  row names: ", row.names(.coords),
					"\n")
	}
	else {
		.coords <- data.frame()
	}
	
	# create instance
	.sml = SensorML(xml = obj, coords = .coords, id = .id, name = .shortName, 
			description = .description, boundedBy = .bb)
	
	if(verbose) cat("[parseSensorML]  Done: ", toString(.sml), "\n")
	
	return(.sml)
}


#
#
#
.smlIdentifier <- function(doc, identifierName, verbose = FALSE) {
	.xpath <- gsub(pattern = .xPathToken, replacement = identifierName,
			x = .smlXPathIdentifier)
	
	if(verbose) cat("[.smlIdentifier] Accessing path ", .xpath, "\n")
	.result <- getNodeSet(doc = doc, path = .xpath,
			namespaces = .sosNamespaceDefinitionsSML)
	
	return(xmlValue(.result[[1]]))
}


#
#
#
plot.SensorML <- function(x, y, ...) {
	.sp <- as(x, "SpatialPointsDataFrame")
	plot(.sp, ...)
}
setMethod("plot", signature(x = "SensorML", y = "missing"),
		plot.SensorML)

