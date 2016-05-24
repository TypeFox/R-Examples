#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


################################################################################
# CLASS:                    REPRESENTATION:
#  setClass                  Classify 'timeSeries'
#  setValidity               Validate 'timeSeries'
#  setMethod                 Initialize 'timeSeries'
################################################################################
# CLASS:                    REPRESENTATION:
#  'signalSeries'            Deprecated S4 Class representation
#  'timeSeries'              Deprecated S4 Class representation
################################################################################


# YC: Note if slots are added or removed, don't forget to edit
# getDataPart,timeSeries-method and setDataPart,timeSeries-method !!


setClass("timeSeries",
   representation(.Data = "matrix",
                  units = "character",
                  positions = "numeric",
                  format = "character",
                  FinCenter = "character",
                  recordIDs = "data.frame",
                  title = "character",
                  documentation = "character"),
   contains = "structure",
   prototype(matrix(NA),
             units = character(0),
             positions = numeric(0),
             format = character(0),
             FinCenter = character(0),
             recordIDs = data.frame(),
             title = character(0),
             documentation = character(0)))


# ------------------------------------------------------------------------------


.validity_timeSeries <- 
function(object) {
    if ((length(object@positions) > 0) &&
        NROW(object) != length(object@positions))
        return("length of '@positions' not equal to '@.Data' extent")
    if (NCOL(object) != length(object@units))
        return("length of '@units' not equal to '@.Data' extent")
    if (NROW(object@recordIDs) > 0 &
        NROW(object@recordIDs) != nrow(object))
        return("length of '@recordIDs' not equal to '@.Data' extent")
  
    # Return Value:
    TRUE
}


setValidity("timeSeries", .validity_timeSeries)


# ------------------------------------------------------------------------------


# Note it is faster to assign manually all slots of the timeSeries objects.
setMethod("initialize", "timeSeries",
    function(.Object,
             .Data = new("matrix"),
             units = character(0),
             positions = numeric(0),
             format = character(0),
             FinCenter = "",
             #<< FIXME: use identical in code rather than FinCenter == ""
             recordIDs = data.frame(),
             title = character(0),
             documentation = character(0))
{

    # as.double -> crucial for speed improvement in subsetting
    if (!is.double(positions)) positions <- as.double(positions)

    .Object <- timeSeries::setDataPart(.Object, value = .Data)
    `slot<-`(.Object, "units", value = units)
    `slot<-`(.Object, "positions", value = positions)
    `slot<-`(.Object, "format", value = format)
    `slot<-`(.Object, "FinCenter", value = FinCenter)
    `slot<-`(.Object, "recordIDs", value = recordIDs)
    `slot<-`(.Object, "title", value = title)
    `slot<-`(.Object, "documentation", value = documentation)

    # Check only one we needs rather than using validObject
    anyStrings <- function(x)
        if (identical(x, TRUE))  character() else x
    error <- anyStrings(.validity_timeSeries(.Object))
    if (length(error) > 0)
        stop(paste("Initialize timeSeries :", error, collapse = "\n"),
             call. = FALSE, domain = NA)

    # Return Value:
    .Object
})


################################################################################


## setClass("signalSeries",
##          representation(
##                         .Data = "matrix",
##                         units = "character",
##                         recordIDs = "data.frame",
##                         title = "character",
##                         documentation = "character"),
##          contains = "structure",
##          validity = function(object) {
##              if (NCOL(getDataPart(object)) != length(object@units))
##                  return("length of '@units' not equal to '@.Data' extent")
##              TRUE
##          })


## # ------------------------------------------------------------------------------


## setClass("timeSeries",
##          representation(positions = "numeric",
##                         format = "character",
##                         FinCenter = "character"),
##          contains = "signalSeries",
##          validity = function(object) {
##              if (NROW(getDataPart(object)) != length(object@positions))
##                  return("length of '@positions' not equal to '@.Data' extent")
##              if (NCOL(getDataPart(object)) != length(object@units))
##                  return("length of '@units' not equal to '@.Data' extent")
##              TRUE
##          })

################################################################################

