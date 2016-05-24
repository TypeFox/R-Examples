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
# Created: 2010-09-17                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#
# samplingTime is the only time that's really used, so set it as default
#
TM_After <- function(propertyName = sosDefaultTempOpPropertyName, time) {
	new("TM_After", propertyName = propertyName, time = time)
}
TM_Before <- function(propertyName = sosDefaultTempOpPropertyName, time) {
	new("TM_Before", propertyName = propertyName, time = time)
}
TM_During <- function(propertyName = sosDefaultTempOpPropertyName, time) {
	new("TM_During", propertyName = propertyName, time = time)
}
TM_Equals <- function(propertyName = sosDefaultTempOpPropertyName, time) {
	new("TM_Equals", propertyName = propertyName, time = time)
}

#
#
#
OgcBBOX <- function(propertyName = sosDefaultSpatialOpPropertyName,
		envelope) {
	new("OgcBBOX", propertyName = propertyName, envelope = envelope)	
}
OgcContains <- function(propertyName = sosDefaultSpatialOpPropertyName,
		geometry = NULL, envelope = NULL) {
	new("OgcContains", propertyName = propertyName, geometry = geometry,
			envelope = envelope)
}
OgcIntersects <- function(propertyName = sosDefaultSpatialOpPropertyName,
		geometry = NULL, envelope = NULL) {
	new("OgcIntersects", propertyName = propertyName, geometry = geometry,
			envelope = envelope)
}
OgcOverlaps <- function(propertyName = sosDefaultSpatialOpPropertyName,
		geometry = NULL, envelope = NULL) {
	new("OgcOverlaps", propertyName = propertyName, geometry = geometry,
			envelope = envelope)
}


#
#
#
setMethod(f = "encodeXML",
		signature = signature(obj = "TM_After", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] TM_After with", toString(obj@time))
			
			.encoded <- .encodeTM(nodeName = ogcTempOpTMAfterName,
					propertyName = obj@propertyName, time = obj@time,
					sos = sos, verbose = verbose)
			return(.encoded)
		}
)
setMethod(f = "encodeXML",
		signature = signature(obj = "TM_Before", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] TM_After with", toString(obj@time),
						"\n")
			
			.encoded <- .encodeTM(nodeName = ogcTempOpTMBeforeName,
					propertyName = obj@propertyName, time = obj@time,
					sos = sos, verbose = verbose)
			return(.encoded)
		}
)
setMethod(f = "encodeXML",
		signature = signature(obj = "TM_During", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] TM_During with", toString(obj@time),
						"\n")
			
			.encoded <- .encodeTM(nodeName = ogcTempOpTMDuringName,
					propertyName = obj@propertyName, time = obj@time,
					sos = sos, verbose = verbose)
			return(.encoded)
		}
)
setMethod(f = "encodeXML",
		signature = signature(obj = "TM_Equals", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] TM_Equals with", toString(obj@time),
						"\n")
			
			.encoded <- .encodeTM(nodeName = ogcTempOpTMEqualsName,
					propertyName = obj@propertyName, time = obj@time,
					sos = sos, verbose = verbose)
			return(.encoded)
		}
)

.encodeTM <- function(nodeName, propertyName, time, sos, verbose = FALSE) {
	if(verbose) cat("[.encodeTM] ", nodeName, "\n")
	
	.tm <- xmlNode(name = nodeName, namespace = ogcNamespacePrefix)
	.pn <- xmlNode(name = ogcPropertyNameName, namespace = ogcNamespacePrefix)
	xmlValue(.pn) <- propertyName
	.tm$children[[1]] <- .pn
	.time <- encodeXML(obj = time, sos = sos, verbose = verbose)
	.tm$children[[2]] <- .time
	
	return(.tm)
}

setMethod(f = "encodeXML",
		signature = signature(obj = "OgcBBOX", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] OgcBBOX with", toString(obj@time),
						"\n")
			
			.bbox <- xmlNode(name = ogcBBOXName, namespace = ogcNamespacePrefix)
			
			.pN <- .createPropertyName(node = .bbox,
					propertyName = obj@propertyName)
			.env <- encodeXML(obj = obj@envelope, sos = sos)
			.bbox <- addChildren(node = .bbox, kids = list(.pN, .env))
			
			return(.bbox)
		}
)

setMethod(f = "encodeXML",
		signature = signature(obj = "OgcContains", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] OgcContains with", toString(obj@time),
						"\n")
			
			.contains <- .encodeBinarySpatialOp(opName = ogcContainsName,
					propertyName = obj@propertyName, geometry = obj@geometry,
					envelope = obj@envelope)
			
			return(.contains)
		}
)

setMethod(f = "encodeXML",
		signature = signature(obj = "OgcIntersects", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] OgcIntersects with",
						toString(obj@time), "\n")
			
			.intersects <- .encodeBinarySpatialOp(opName = ogcIntersectsName,
					propertyName = obj@propertyName, geometry = obj@geometry,
					envelope = obj@envelope)
			
			return(.intersects)
		}
)

setMethod(f = "encodeXML",
		signature = signature(obj = "OgcOverlaps", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] OgcOverlaps with", toString(obj@time),
						"\n")
			
			.overlaps <- .encodeBinarySpatialOp(opName = ogcOverlapsName,
					propertyName = obj@propertyName, geometry = obj@geometry,
					envelope = obj@envelope, sos = sos)
			
			return(.overlaps)
		}
)

.encodeBinarySpatialOp <- function(opName, propertyName, geometry, envelope,
		sos) {
	.spOp <- xmlNode(name = opName, namespace = ogcNamespacePrefix)
	
	.pN <- .createPropertyName(node = .spOp,
			propertyName = propertyName)
	
	# switch between geometry and envelope
	if(!is.null(geometry)) {
		.geomOrEnv <- encodeXML(obj = geometry, sos = sos)
	}
	else if(!is.null(envelope)) {
		.geomOrEnv <- encodeXML(obj = envelope, sos = sos)
	}
	else {
		warning("At least one of geometry or envelope has to be set.")
		.geomOrEnv <- NULL
	}
	
	.spOp <- addChildren(node = .spOp, kids = list(.pN, .geomOrEnv))
	
	return(.spOp)
}

.createPropertyName <- function(node, propertyName) {
	.pN <- xmlNode(name = ogcPropertyNameName, namespace = ogcNamespacePrefix)
	xmlValue(.pN) <- propertyName
	return(.pN)
}

setMethod(f = "encodeXML",
		signature = signature(obj = "OgcComparisonOps", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeXML] OgcComparisonOps with",
						toString(obj@time), "\n")
			warning("Encoding of OgcComparisonOps not implemented yet! Returning obj as is...")
			return(obj)
		}
)

#
# see: http://www.oostethys.org/best-practices/best-practices-get
#
setMethod(f = "encodeKVP",
		signature = signature(obj = "OgcBinaryTemporalOp", sos = "SOS"),
		def = function(obj, sos, verbose) {
			if(verbose) cat("[encodeKVP] temporalOps: ", toString(obj), "\n")
			.time <- NULL
			.tempOpTime <- obj@time
			
			if(class(.tempOpTime) == "GmlTimeInstant") {
				if(verbose)
					cat("[encodeKVP] Encoding instant.\n")
				.time <- encodeKVP(.tempOpTime@timePosition@time, sos = sos,
						verbose = verbose)
			}
			# ignore type, because temporal operators are not supportded by the
			# GET binding
			else if (class(.tempOpTime) == "GmlTimePeriod") {
				if(!is.null(.tempOpTime@begin) && !is.null(.tempOpTime@end)) {
					if(verbose)
						cat("[encodeKVP] Encoding period with begin and end.\n")
					.begin <- encodeKVP(.tempOpTime@begin@time@timePosition, sos = sos,
							verbose = verbose)
					.end <- encodeKVP(.tempOpTime@end@time@timePosition, sos = sos,
							verbose = verbose)
					.time <- paste(.begin, "/", .end, sep = "")
				}
				else if(!is.null(.tempOpTime@beginPosition)
						&& !is.null(.tempOpTime@endPosition)) {
					if(verbose)
						cat("[encodeKVP] Encoding period with beginPosition and endPosition.\n")
					.begin <- encodeKVP(.tempOpTime@beginPosition@time, sos = sos,
							verbose = verbose)
					.end <- encodeKVP(.tempOpTime@endPosition@time, sos = sos,
							verbose = verbose)
					.time <- paste(.begin, "/", .end, sep = "")
				}
				else {
					stop(paste("Incomplete gml:TimePeriod:",
									toString(.tempOpTime)))
				}
			}
			else {
				stop(paste("Cannot encode given object as KVP",
								toString(.tempOpTime)))
			}
			
			return(.time)
		}
)
