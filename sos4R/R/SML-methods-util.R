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
# Created: 2011-02-11                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################


#
#
#
setMethod(f = "sosId", signature = signature(obj = "SensorML"),
		def = function(obj) {
			return(obj@id)
		})

#
#
#
setMethod(f = "sosName", signature = signature(obj = "SensorML"),
		def = function(obj) {
			return(obj@name)
		})

#
#
#
setMethod(f = "sosAbstract", signature = signature(obj = "SensorML"),
		def = function(obj) {
			return(obj@description)
		})

#
#
#
setMethod(f = "sosGetCRS",
		signature = c(obj = "SensorML"),
		def = function(obj, verbose = FALSE) {
			.coords <- sosCoordinates(obj, handleNames = FALSE)
			.crs <- sosGetCRS(attributes(.coords)[["referenceFrame"]],
					verbose = verbose)
			return(.crs)
		}
)

#
#
#
setMethod(f = "sosBoundedBy",
		signature = signature(obj = "SensorML"),
		def = function(obj, sos, verbose = FALSE) {
			return(obj@boundedBy)
		})


#
# extract the coordinates from the SensorML document and return as a data.frame
#
setMethod(f = "sosCoordinates", signature = signature(obj = "SensorML"),
		def = function(obj, handleNames = TRUE, verbose = FALSE) {
			.coords <- obj@coords
			
			if(handleNames) {
				if(verbose) cat("Handling coordinate names enabled!\n")
				
				# handle coord name variations
				.coordsNames <- names(.coords)
				
				# easting/northing
				if(any(.coordsNames == "easting")) {
					.coordsNames[which(.coordsNames == "easting")] <- "x"
					warning("Changed coordinate name for 'easting' to 'x'.")
				}
				if(any(.coordsNames == "northing")) {
					.coordsNames[which(.coordsNames == "northing")] <- "y"
					warning("Changed coordinate name for 'northing' to 'y'.")
				}
				
				# longitude/latitude
				if(any(.coordsNames == "longitude")) {
					.coordsNames[which(.coordsNames == "longitude")] <- "x"
					warning("Changed coordinate name for 'longitude' to 'x'.")
				}
				if(any(.coordsNames == "latitude")) {
					.coordsNames[which(.coordsNames == "latitude")] <- "y"
					warning("Changed coordinate name for 'latitude' to 'y'.")
				}
				
				# elevation/altitude
				if(any(.coordsNames == "elevation")) {
					.coordsNames[which(.coordsNames == "elevation")] <- "z"
					warning("Changed coordinate name for 'elevation' to 'z'.")
				}
				if(any(.coordsNames == "altitude")) {
					.coordsNames[which(.coordsNames == "altitude")] <- "z"
					warning("Changed coordinate name for 'altitude' to 'z'.")
				}
				
				names(.coords) <- .coordsNames
				#cat(.coordsNames)
			}
			
			return(.coords)
		}
)
