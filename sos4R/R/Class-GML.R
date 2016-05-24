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
# Created: 2010-09-08                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#
# Classes are based on GML
# 
# http://www.opengeospatial.org/standards/gml
#

#
# abstract super classes of the actually needed TimeInstant and TimePeriod
#
setClass("GmlTimeObject",
		representation(# optional:
				id = "character"),
		contains = c("VIRTUAL"),
		validity = function(object) {
			#print("Entering validation: GmlTime")
			# TODO implement validity function
			return(TRUE)
		}
)
setClassUnion(name = "GmlTimeObjectOrNULL",
		members = c("GmlTimeObject", "NULL"))

setClass("GmlTimePrimitive",
		representation(# optional:
				relatedTimes = "list"),
		prototype = list(relatedTimes = list()),
		contains = c("VIRTUAL", "GmlTimeObject"),
		validity = function(object) {
			#print("Entering validation: GmlTimePrimitive")
			# TODO implement validity function
			# related times must be GmlTimePrimitive
			return(TRUE)
		}
)

setClass("GmlTimeGeometricPrimitive",
		representation(# optional:
				frame = "character"),
		contains = c("VIRTUAL", "GmlTimePrimitive"),
		validity = function(object) {
			#print("Entering validation: GmlTimeGeometricPrimitive")
			# TODO implement validity function
			return(TRUE)
		}
)

#
# GmlTimePosition wraps a POSIXct
#
setClass("GmlTimePosition",
		representation(time = "POSIXt",
				# optional:
				frame = "character",
				calendarEraName = "character",
				indeterminatePosition = "character"
		),
		prototype = list(time = as.POSIXct(NA)),
		validity = function(object) {
			#print("Entering validation: GmlTimePosition")
			# TODO implement validity function
			# time needs to be set			
			return(TRUE)
		}		
)
setClassUnion(name = "GmlTimePositionOrNULL",
		members = c("GmlTimePosition", "NULL"))

#
#
#
setClass("GmlTimeInstant",
		representation(timePosition = "GmlTimePosition"),
		prototype = list(timePosition = NULL),
		contains = "GmlTimeGeometricPrimitive",
		validity = function(object) {
			#print("Entering validation: GmlTimeInstant")
			# TODO implement validity function
			# timePosition needs to be set			
			return(TRUE)
		}		
)
setClassUnion(name = "GmlTimeInstantOrNULL",
		members = c("GmlTimeInstant", "NULL"))

#
#
#
setClass("GmlTimeInstantProperty",
		representation(href = "character", time = "GmlTimeInstantOrNULL"),
		prototype = list(href = as.character(NA), time = NULL),
		contains = "GmlTimeGeometricPrimitive",
		validity = function(object) {
			#print("Entering validation: GmlTimeInstant")
			# TODO implement validity function
			# time needs to be set			
			return(TRUE)
		}
)
setClassUnion(name = "GmlTimeInstantPropertyOrNULL",
		members = c("GmlTimeInstantProperty", "NULL"))

#
#
#
setClass("GmlTimeInterval",
		representation(interval = "character",
				unit = "character",
				# optional:
				radix = "integer",
				factor = "integer"
		),
		prototype = list(interval = as.character(NA), unit = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: GmlTimeInterval")
			# TODO implement validity function
			# interval and unit need to be set, radix must be positive
			return(TRUE)
		}		
)
setClassUnion(name = "GmlTimeIntervalOrNULL",
		members = c("GmlTimeInterval", "NULL"))

#
#
#
setClass("GmlTimePeriod",
		representation(begin = "GmlTimeInstantPropertyOrNULL",
				beginPosition = "GmlTimePositionOrNULL",
				end = "GmlTimeInstantPropertyOrNULL",
				endPosition = "GmlTimePositionOrNULL",
				# optional:
				# for brevity, the TimeDurationType layer is removed here
				duration = "character",
				timeInterval = "GmlTimeIntervalOrNULL"
		),
		prototype = list(begin = NULL, beginPosition = NULL, end = NULL, 
				endPosition = NULL),
		contains = "GmlTimeGeometricPrimitive",
		validity = function(object) {
			#print("Entering validation: GmlTimeInstant")
			# TODO implement validity function
			# either both begin and end, or beginPosition and endPosition need to be set.
			# only one of the optional duration and timeInterval can be set!
			return(TRUE)
		}		
)

#
#
#
setClass("GmlFeature",
		representation(id = "character"),
		contains = "VIRTUAL",
		#prototype = list(),
		validity = function(object) {
			#print("Entering validation: GmlFeature")
			return(TRUE)
		}
)
setClassUnion(name = "GmlFeatureOrNULL", members = c("GmlFeature", "NULL"))

#
#
#
setClass("GmlFeatureProperty",
		representation(href = "character",	
				feature = "GmlFeatureOrNULL"),
		prototype = list(href = as.character(NA), feature = NULL),
		validity = function(object) {
			#print("Entering validation: GmlFeatureProperty")
			# TODO implement validity function
			# one of parameters has to be set
			return(TRUE)
		}
)
setClassUnion(name = "GmlFeatureOrGmlFeaturePropertyOrNULL",
	members = c("GmlFeatureProperty", "GmlFeature", "NULL"))


#
#
#
setClass("GmlFeatureCollection",
		representation(featureMembers = "list"),
		contains = c("GmlFeature"),
		validity = function(object) {
			#print("Entering validation: GmlFeatureCollection")
			return(TRUE)
		}
)

#
#
#
setClass("GmlDirectPosition",
		representation(pos = "character",
				# optional:
				srsName = "character",
				srsDimension = "integer",
				axisLabels = "character",
				uomLabels = "character"),
		prototype = list(pos = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: GmlDirectPosition")
			# TODO implement validity function
			# pos string is required
			return(TRUE)
		}
)
setClassUnion(name = "GmlDirectPositionOrNULL",
		members = c("GmlDirectPosition", "NULL"))

#
#
#
setClass("GmlGeometry",
		representation(id = "character",
				srsName = "character",
				srsDimension = "integer",
				axisLabels = "character",
				uomLabels = "character"),
		contains = c("VIRTUAL"),
		validity = function(object) {
			#print("Entering validation: GmlGeometry")
			# TODO implement validity function
			# one of parameters has to be set
			return(TRUE)
		}
)



#
# gml:coordinates and gml:coord are deprecated
#
setClass("GmlPoint",
		representation(pos = "GmlDirectPosition"),
		prototype = list(pos = NULL),
		contains = c("GmlGeometry"),
		validity = function(object) {
			#print("Entering validation: GmlPoint")
			# TODO implement validity function
			# pos must be set
			return(TRUE)
		}
)
setClassUnion(name = "GmlPointOrNULL", members = c("GmlPoint", "NULL"))

#
# NOT IMPLEMENTED
#
setClass("GmlLineString",
		representation(poss = "list", points = "list", posList = "ANY"),
		prototype = list(),
		contains = c("GmlGeometry"),
		validity = function(object) {
			#print("Entering validation: GmlLineString")
			# TODO implement validity function
			# one of the possible lists or posList has to be set
			return(TRUE)
		}
)

#
# NOT IMPLEMENTED
#
setClass("GmlPolygon",
		representation(exterior = "ANY", interior = "ANY"),
		prototype = list(),
		contains = c("GmlGeometry"),
		validity = function(object) {
			#print("Entering validation: GmlPolygon")
			# TODO implement validity function
			return(TRUE)
		}
)



#
#
#
setClass("GmlPointProperty",
		representation(href = "character",	
				point = "GmlPointOrNULL"),
		#prototype = list(),
		validity = function(object) {
			#print("Entering validation: GmlPointProperty")
			# TODO implement validity function
			# one of parameters has to be set, point has to be NA or GmlPoint
			return(TRUE)
		}
)

#
#
#
setClass("GmlEnvelope",
		representation(
				lowerCorner = "GmlDirectPosition",
				upperCorner = "GmlDirectPosition",
				srsName = "character",
				srsDimension = "integer",
				axisLabels = "character",
				uomLabels = "character"),
		prototype = list(lowerCorner = NULL, upperCorner = NULL),
		validity = function(object) {
			#print("Entering validation: GmlEnvelope")
			# TODO implement validity function
			# srsDimension must be positive
			return(TRUE)
		}
)

################################################################################
#
#
#
#
setClass("GmlMeasure",
		representation(value = "numeric",
				uom = "character"),
		prototype = list(value = as.numeric(NA), uom = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: GmlMeasure")
			# TODO implement validity function
			# both parameters are required!
			return(TRUE)
		}
)
