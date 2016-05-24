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
# Classes are based on Observations & Measurements - Part 2 - Sampling Features
# 
# http://www.opengeospatial.org/standards/om
#

#
# 
#
setClass("SaSamplingPoint",
		representation(sampledFeatures = "list",
				position = "GmlPointProperty",
				# optional:
				relatedObservation = "list",
				relatedSamplingFeature = "list",
				surveyDetails = "ANY"),
		prototype = list(sampledFeatures = list(NA), position = NULL),
		contains = "GmlFeature",
		validity = function(object) {
			#print("Entering validation: SaSamplingPoint")
			# TODO implement validity function
			# sampledFeatures list must contain > 0 gml:_Feature instances
			# related observations must be OmObservationProperty
			return(TRUE)
		}
)

#
#
#
setClass("SaSamplingSurface",
		representation(sampledFeatures = "list",
				shape = "ANY",
				# optional:
				relatedObservation = "list",
				relatedSamplingFeature = "list",
				surveyDetails = "ANY",
				area = "ANY"),
		prototype = list(sampledFeatures = list(NA), shape = NULL),
		contains = "GmlFeature",
		validity = function(object) {
			#print("Entering validation: SaSamplingSurface")
			# TODO implement validity function
			# sampledFeatures list must contain > 0 gml:_Feature instances
			# related observations must be of type OmObservationProperty
			return(TRUE)
		}
)

