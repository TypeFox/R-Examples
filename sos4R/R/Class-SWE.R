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
setClass("SwePhenomenon",
		representation(id = "character", name = "character",
				# optional:
				description = "character"),
		prototype = list(id = as.character(NA), name = as.character(NA)),
		validity = function(object) {
			#print("Entering validation: SwePhenomenon")
			# TODO implement validity function
			# one of parameters has to be set
			return(TRUE)
		}
)
setClassUnion(name = "SwePhenomenonOrNULL",
		members = c("SwePhenomenon", "NULL"))

#
#
#
setClass("SwePhenomenonProperty",
		representation(href = "character", phenomenon = "SwePhenomenonOrNULL"),
		prototype = list(href = as.character(NA), phenomenon = NULL),
		validity = function(object) {
			#print("Entering validation: SwePhenomenonProperty")
			# TODO implement validity function
			# one of parameters has to be set, phenomenon has to be SwePhenomenon if set
			return(TRUE)
		}
)
setClassUnion(name = "SwePhenomenonPropertyOrNULL",
		members = c("SwePhenomenonProperty", "NULL"))

#
# 52N SOS only supports/returns CompositePhenomenon, the intermediate
# CompoundPhenomenon is intentionally left out for brevity
#
setClass("SweCompositePhenomenon",
		representation(dimension = "integer", components = "list",
				# optional: 
				base = "SwePhenomenonPropertyOrNULL"),
		prototype = list(dimension = NA_integer_, components = list(NA)),
		contains = "SwePhenomenon",
		validity = function(object) {
			#print("Entering validation: SweCompositePhenomenon")
			# TODO implement validity function
			# components needs to be a list of SwePhenomenonProperty instances
			return(TRUE)
		}
)

#
#
#
setClass("SweTextBlock",
		representation(tokenSeparator = "character",
				blockSeparator = "character",
				decimalSeparator = "character",
				#optional:
				id = "character"),
		prototype = list(tokenSeparator = NA_character_,
				blockSeparator = NA_character_,
				decimalSeparator = NA_character_),
		validity = function(object) {
			#print("Entering validation: SweTextBlock")
			# TODO implement validity function
			return(TRUE)
		}
)


################################################################################
#
# Other SOS related SWE elements, like DataRecord, SimpleDataRecord, values, and
# so forth are parsed directly into R structures.
#

