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
# Created: 2010-10-19                                                          #
# Project: sos4R - visit the project web page, http://www.nordholmen.net/sos4r #
#                                                                              #
################################################################################

#
#
#
as.data.frame.OmObservation = function(x, row.names, optional, ...) {
	return(sosResult(x))
}
setAs(from = "OmObservation", to = "data.frame",
		def = function(from) {
			as.data.frame.OmObservation(from)
		}
)

#
#
#
as.data.frame.OmMeasurement = function(x, row.names, optional, ...) {
	return(sosResult(x))
}
setAs(from = "OmMeasurement", to = "data.frame",
		def = function(from) {
			as.data.frame.OmMeasurement(from)
		}
)

#
#
#
as.list.OmObservationCollection = function(x, ...) {
	return(x@members)
}
setAs(from = "OmObservationCollection", to = "list",
		def = function(from) {
			as.list.OmObservationCollection(from)
		}
)

#
#
#
as.SpatialPointsDataFrame.OmObservationCollection = function(x, ...) {
	.result <- sosResult(x, coordinates = TRUE)
	.crs <- sosGetCRS(x)
	
	if(length(.crs) > 1)
		stop("Spatial Reference System is not unambiguous, cannot convert.")
	
	.spdf <- .resultDataFrameToSpatialPointsDataFrame(result = .result,
			crs = .crs)
	return(.spdf)
}
setAs(from = "OmObservationCollection", to = "SpatialPointsDataFrame",
		def = function(from) {
			as.SpatialPointsDataFrame.OmObservationCollection(from)
		}
)

#
#
#
as.SpatialPointsDataFrame.OmObservation = function(x, ...) {
	.crs <- sosGetCRS(x)
	.result <- sosResult(x, coordinates = TRUE)
	
	.spdf <- .resultDataFrameToSpatialPointsDataFrame(result = .result,
			crs = .crs)
	return(.spdf)
}
setAs(from = "OmObservation", to = "SpatialPointsDataFrame",
		def = function(from) {
			as.SpatialPointsDataFrame.OmObservation(from)
		}
)

#
#
#
as.SpatialPointsDataFrame.OmMeasurement = function(x, ...) {
	.crs <- sosGetCRS(x)
	.result <- sosResult(x, coordinates = TRUE)
	
	.spdf <- .resultDataFrameToSpatialPointsDataFrame(result = .result,
			crs = .crs)
	return(.spdf)
}
setAs(from = "OmObservation", to = "SpatialPointsDataFrame",
		def = function(from) {
			as.SpatialPointsDataFrame.OmObservation(from)
		}
)

################################################################################
#
.resultDataFrameToSpatialPointsDataFrame <- function(result, crs) {
	# TODO fix column order, which is x~y according to ?coordinates
	.coordCols <- match(c(sosDefaultColumnNameLat, sosDefaultColumnNameLon),
			colnames(result))
	
	.spdf <- SpatialPointsDataFrame(
			coords = result[, .coordCols],
			data = result[, -.coordCols],
			proj4string = crs)
	
	return(.spdf)
}

