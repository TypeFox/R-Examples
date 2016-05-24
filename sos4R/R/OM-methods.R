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
# construction methods
#
OmObservationCollection <- function(members, boundedBy) {
	new("OmObservationCollection", members = members, boundedBy = boundedBy)
}

OmObservation <- function(samplingTime, procedure, observedProperty,
		featureOfInterest, result, metadata = NA, resultTime = NULL,
		resultQuality = NA,	parameter = NA) {
	new("OmObservation", samplingTime = samplingTime, procedure = procedure,
			observedProperty = observedProperty,
			featureOfInterest = featureOfInterest, result = result,
			metadata = metadata, resultTime = resultTime,
			resultQuality = resultQuality,
			parameter = parameter)
}

OmObservationProperty <- function(href = as.character(NA), obs = NULL) {
	new("OmObservationProperty", href = href, obs = obs)
}

OmMeasurement <- function(samplingTime, procedure, observedProperty,
		featureOfInterest, result, metadata = NA, resultTime = NULL,
		resultQuality = NA,	parameter = NA) {
	new("OmMeasurement", samplingTime = samplingTime, procedure = procedure,
			observedProperty = observedProperty,
			featureOfInterest = featureOfInterest, result = result,
			metadata = metadata, resultTime = resultTime,
			resultQuality = resultQuality,
			parameter = parameter)
}


################################################################################
#
# Some problem with this function: Could not find function "getGeneric" ...
#setMethod(f = "length", signature = signature(x = "OmObservationCollection"),
#		def = function(x) {
#			.l <- length(x@members)
#			return(.l)
#		}
#)
length.OmObservationCollection <- function(x) {
	length(x@members)
}

setMethod(f = "[[", signature = signature(x = "OmObservationCollection",
				i = "ANY", j = "missing"), 
		def = function(x, i, j, ...) {
			if(is.numeric(i)) {
				return(x@members[[i]])
			}
			else {
				warning("Indexing only supported with numeric values!")
			}
		}
)

.getObservationsWithObservedProperty <- function(coll, obsProp) {
	.obsProperties <- sosObservedProperties(coll)

	if(any(is.na(.obsProperties))) {
		# warning("NA values in observed property list.")
		# remove NAs
		.obsProperties <- .obsProperties[which(!is.na(.obsProperties))]
	}
	
	if(length(.obsProperties) < 1)
		return(list())
	
	.idx <- c()
	
	for (i in seq(1:length(.obsProperties))) {
		if(is.list(.obsProperties[[i]])) {
			.current <- .obsProperties[[i]]
#			cat(i, ": current:", .current, "\n")
			if(any(.current == obsProp)) {
				.idx <- c(.idx, i)
#				cat("found index: ", i, ": ", .idx, "\n")
			}
		}
		else {
			if(.obsProperties[[i]] == obsProp) {
				.idx <- c(.idx, i)
#				cat("found index: ", i, ": ", .idx, "\n")
			}
		}
	}
#	cat("Found observed property ", obsProp, " at indices", .idx, "\n")
	if(length(.idx) == 0)
		return(list())
	else
		return(coll[.idx])
}
.getObservationsWithProcedure <- function(coll, procedure) {
	.procedures <- sosProcedures(coll)
	.idx <- which(.procedures %in% procedure)
#	cat("Found procedure ", procedure, " at indices", .idx, "\n")
	if(length(.idx) == 0)
		return(list())
	else
		return(coll[.idx])
}
.getObservationsWithFoiId <- function(coll, foiId) {
	.featureIds <- sosFeatureIds(coll)
	.idx <- which(.featureIds %in% foiId)
#	cat("Found foi ", foiId, " at indices", .idx, "\n")
	if(length(.idx) == 0)
		return(list())
	else
		return(coll[.idx])
}

setMethod(f = "[", signature = signature(x= "OmObservationCollection", 
				i = "ANY", j = "ANY"),
		def = function(x, i, j, ...) {
			if (missing(j)) {
				if(is.numeric(i)) {
					return(x@members[i])
				}
				else {
#					cat("Try subsetting with", i, "\n")
					# subset the collection by procedure or observed property
					.byProc <- .getObservationsWithProcedure(x, i)
#					cat("by procedures: ", toString(.byProc), "\n")
					if(length(.byProc) > 0)
						return(.byProc)
					
					.byObsProp <- .getObservationsWithObservedProperty(x, i)
#					cat("by obs prop: ", toString(.byObsProp), "\n")
					if(length(.byObsProp) > 0)
						return(.byObsProp)
					
					.byFoiId <- .getObservationsWithFoiId(x, i)
#					cat("by foi id: ", toString(.byObsProp), "\n")
					if(length(.byFoiId) > 0)
						return(.byFoiId)
					
					return(list())
				}
			}
			else return(x@members[i:j])
		}
)

################################################################################
#
names.OmObservation <- function(x) {
	.name <- paste(sosProcedures(x), sosObservedProperties(x), sosFeatureIds(x),
			sep = "_")
	names(.name) <- "proc_obsProp_foiID"
	return(.name)
}

names.OmMeasurement <- function(x) {
	.name <- paste(sosProcedures(x), sosObservedProperties(x), sosFeatureIds(x),
			sep = "_")
	names(.name) <- "proc_obsProp_foiID"
	return(.name)
}

names.OmObservationCollection <- function(x) {
	.names <- sapply(x@members, names)
	return(.names)
}
