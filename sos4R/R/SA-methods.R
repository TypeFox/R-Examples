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
# Created: 2010-09-09                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#
# construction methods
#
SaSamplingPoint <- function(sampledFeatures, position,
		relatedObservation = list(NA), relatedSamplingFeature = list(NA),
		surveyDetails = NA, id = NA_character_) {
	new("SaSamplingPoint", sampledFeatures = sampledFeatures,
			position = position, relatedObservation = relatedObservation,
			surveyDetails = surveyDetails, id = id)
}

#
#
#
SaSamplingSurface <- function(sampledFeatures, shape,
		relatedObservation = list(NA), relatedSamplingFeature = list(NA),
		surveyDetails = NA, position = NA) {
	new("SaSamplingSurface", sampledFeatures = sampledFeatures,
			shape = shape, relatedObservation = relatedObservation,
			relatedSamplingFeature = relatedSamplingFeature, 
			surveyDetails = surveyDetails, position = position)
}
