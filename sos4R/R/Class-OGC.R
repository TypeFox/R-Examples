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
# TEMPORAL

#
# ogc:temporalOps needed in GetObservation Request:
# 52N SOS only supports TimeInstant and TimePeriod for eventTime, see
# HttpPostRequestDecoder.java, so time can be a GmlTimeGeometricPrimitive
#
setClass("OgcBinaryTemporalOp",
		representation(propertyName = "character",
				time = "GmlTimeGeometricPrimitive"),
		prototype = list(propertyName = as.character(NA), time = NULL),
		contains = c("VIRTUAL"),
		validity = function(object) {
			#print("Entering validation: OgcBinaryTemporalOp")
			# TODO implement validity function
			return(TRUE)
		}
)
setClassUnion(name = "OgcBinaryTemporalOpOrNULL",
		members = c("OgcBinaryTemporalOp", "NULL"))

#
# after and before are only allowed with time instant
#
setClass("TM_After",
		representation(time = "GmlTimeInstant"),
		contains = c("OgcBinaryTemporalOp"),
		validity = function(object) {
			#print("Entering validation: TM_After")
			# TODO implement validity function
			return(TRUE)
		}
)
setClass("TM_Before",
		representation(time = "GmlTimeInstant"),
		contains = c("OgcBinaryTemporalOp"),
		validity = function(object) {
			#print("Entering validation: TM_Before")
			# TODO implement validity function
			return(TRUE)
		}
)

#
# during makes only sense with time period
#
setClass("TM_During",
		representation(time = "GmlTimePeriod"),
		contains = c("OgcBinaryTemporalOp"),
		validity = function(object) {
			#print("Entering validation: TM_During")
			# TODO implement validity function
			return(TRUE)
		}
)

#
# equals allows both time instant and period
#
setClass("TM_Equals",
		contains = c("OgcBinaryTemporalOp"),
		validity = function(object) {
			#print("Entering validation: TM_Equals")
			# TODO implement validity function
			return(TRUE)
		}
)


################################################################################
# SPATIAL

#
#
#
setClass("OgcSpatialOps",
		contains = c("VIRTUAL"),
		validity = function(object) {
			#print("Entering validation: OgcSpatialOps")
			return(TRUE)
		}
)
setClassUnion(name = "OgcSpatialOpsOrNULL",
		members = c("OgcSpatialOps", "NULL"))

setClass("OgcBBOX",
		representation(propertyName = "character",
				envelope = "GmlEnvelope"),
		contains = c("OgcSpatialOps"),
		validity = function(object) {
			#print("Entering validation: OgcBBOX")
			return(TRUE)
		}
)

#
#
#
setClass("OgcBinarySpatialOp",
		representation(propertyName = "character",
				geometry = "GmlGeometry",
				envelope = "GmlEnvelope"),
		contains = c("VIRTUAL", "OgcSpatialOps"),
		prototype = list(propertyName = as.character(NA), geometry = NULL,
				envelope = NULL),
		validity = function(object) {
			print("Entering validation: OgcBinarySpatialOp")
			# TODO implement validity function
			# only one of geometry of envelope can be set
			return(TRUE)
		}
)
setClass("OgcContains",
		contains = c("OgcBinarySpatialOp"),
		validity = function(object) {
			print("Entering validation: OgcContains")
			return(TRUE)
		}
)
setClass("OgcIntersects",
		contains = c("OgcBinarySpatialOp"),
		validity = function(object) {
			print("Entering validation: OgcIntersects")
			return(TRUE)
		}
)
setClass("OgcOverlaps",
		contains = c("OgcBinarySpatialOp"),
		validity = function(object) {
			print("Entering validation: OgcOverlaps")
			return(TRUE)
		}
)


################################################################################
# RESULT FILTERING

#
#
#
setClass("OgcComparisonOps",
		#contains = c("VIRTUAL"),
		validity = function(object) {
			#print("Entering validation: OgcSpatialOps")
			return(TRUE)
		}
)
#setClassUnion(name = "OgcComparisonOpsOrXMLOrNULL",
#		members = c("OgcComparisonOps", "XMLNode", "NULL", "XMLAbstractNode"))
#"XMLPINode", "XMLCommentNode", "XMLProcessingInstruction",
#"XMLCDataNode", "RXMLAbstractNode", "XMLHashTreeNode", "XMLTextNode",
#"XMLPINode", "XMLCommentNode", "XMLProcessingInstruction",
#"XMLCDataNode", "XMLAttributeNode"
#
# Removed this class union to avoid warnings on installation:
# - DONE: manually check in validity function, not so nice: https://stat.ethz.ch/pipermail/bioc-devel/2010-August/002292.html
# - follow up on this old thread: http://www.mail-archive.com/r-devel@r-project.org/msg15088.html
# - another old thread that went unanswered: http://tolstoy.newcastle.edu.au/R/e2/devel/06/12/1328.html
