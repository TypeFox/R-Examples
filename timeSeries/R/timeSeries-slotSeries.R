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
# FUNCTION:                        DESCRIPTION:
#  series,timeSeries                Get data slot from a'timeSeries'  
#  series<-,timeSeries,ANY          Set new data slot to a 'timeSeries'  
#  series<-,timeSeries,matrix       Set new data slot to a 'timeSeries'  
# SYNONYMES:                       DESCRIPTION:
#  coredata,timeSeries              Get data slot from a'timeSeries'  
#  coredata<-,timeSeries,ANY        Set new data slot to a 'timeSeries'  
#  coredata<-,timeSeries,matrix     Set new data slot to a 'timeSeries'  
# FUNCTION:                        DESCRIPTION:
#  getSeries
#  getSeries.default
#  getSeries.timeSeries             Get data slot from a 'timeSeries'  
#  setSeries<-                      Set new data slot to a 'timeSeries' 
################################################################################
# DEPRECATED:                      DESCRIPTION:
#  seriesData                       Deprecated, use series
################################################################################


setMethod("series", "timeSeries", 
    function(x)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #    Returns the series Data of an ordered data object.

    # Arguments:
    #   x - a 'timeSeries' object

    # Value:
    #    Returns an object of class 'matrix'.

    # FUNCTION:

    # Get Data Slot:
    ans <- as.matrix(x)

    # Return Value:
    ans
}
)


# ------------------------------------------------------------------------------


setMethod("series<-", signature(x = "timeSeries", value = "ANY"), 
    function(x, value)
{
    # A function implemented by Yohan Chalabi
    
    # Return Value:
    callGeneric(x, as(value, "matrix"))
}
)


# ------------------------------------------------------------------------------


setMethod("series<-", signature(x = "timeSeries", value = "matrix"),
    function(x, value)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #    Assign the series Data to a timeSeries object.

    # Arguments:
    #   object - a 'timeSeries' object

    # Value:
    #    Assign to be assign as series Data of a timeSeries.

    # FUNCTION:

    # if value same dimension as time series
    # we we can assign the value directly to @.Data
    # This can speed up math Ops significantly
    if (identical(dim(x), dim(value))) {
        x@.Data <- value
        if (!is.null(cn <- colnames(value)))
            colnames(x) <- cn
        return(x)
    }

    if (is.null(charvec <- rownames(value)))
        charvec <- rownames(x)
    if (is.null(units <- colnames(value)))
        units <- colnames(value)

    # now that we have charvec and units, better to remove
    # dimnames of value to avoid problems
    attr(value, "dimnames") <- NULL

    if (!identical(length(units), NCOL(value)))
        units <- NULL

    # if now same dim , drop charvec and returns .signalSeries
    if (!identical(length(charvec), NROW(value)))
        return(.signalSeries(value, units))

    format <- x@format
    zone <- FinCenter <- finCenter(x)
    title <- x@title
    documentation <- x@documentation
    recordIDs <-
        if (identical(NROW(x), NROW(value)))
            x@recordIDs
        else
            data.frame()

    # Return Value:
    timeSeries(data = value,
        charvec = charvec,
        units = units,
        format = format,
        zone = zone,
        FinCenter = FinCenter,
        recordIDs = recordIDs,
        title = title)
}
)


################################################################################
# COREDATA SYNONYME


setMethod("coredata", "timeSeries", 
    function(x)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #    Returns the series Data of an ordered data object.

    # Arguments:
    #   x - a 'timeSeries' object

    # Value:
    #    Returns an object of class 'matrix'.

    # FUNCTION:

    # Get Data Slot:
    ans <- as.matrix(x)

    # Return Value:
    ans
})


# ------------------------------------------------------------------------------


setMethod("coredata<-", signature(x = "timeSeries", value = "ANY"), 
    function(x, value)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi
    
    # Return Value:
    callGeneric(x, as(value, "matrix"))
})


# ------------------------------------------------------------------------------


setMethod("coredata<-", signature(x = "timeSeries", value = "matrix"),
    function(x, value)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #    Assign the series Data to a timeSeries object.

    # Arguments:
    #   object - a 'timeSeries' object

    # Value:
    #    Assign to be assign as series Data of a timeSeries.

    # FUNCTION:

    # if value same dimension as time series
    # we we can assign the value directly to @.Data
    # This can speed up math Ops significantly
    if (identical(dim(x), dim(value))) {
        x@.Data <- value
        if (!is.null(cn <- colnames(value)))
            colnames(x) <- cn
        return(x)
    }

    if (is.null(charvec <- rownames(value)))
        charvec <- rownames(x)
    if (is.null(units <- colnames(value)))
        units <- colnames(value)

    # now that we have charvec and units, better to remove
    # dimnames of value to avoid problems
    attr(value, "dimnames") <- NULL

    if (!identical(length(units), NCOL(value)))
        units <- NULL

    # if now same dim , drop charvec and returns .signalSeries
    if (!identical(length(charvec), NROW(value)))
        return(.signalSeries(value, units))

    format <- x@format
    zone <- FinCenter <- finCenter(x)
    title <- x@title
    documentation <- x@documentation
    recordIDs <-
        if (identical(NROW(x), NROW(value)))
            x@recordIDs
        else
            data.frame()

    # Return Value:
    timeSeries(data = value,
        charvec = charvec,
        units = units,
        format = format,
        zone = zone,
        FinCenter = FinCenter,
        recordIDs = recordIDs,
        title = title)
})


################################################################################


## getSeries <-
##     function(x)
## {
##     # A function implemented by Diethelm Wuertz
##     
##     # Note:
##     #   Used for getSeries methods in fPortfolio package.
## 
##     # FUNCTION: 
##     
##     # Return Value:
##     UseMethod("getSeries")
## }


# ------------------------------------------------------------------------------


## getSeries.default <- 
##     function(x)
## {
##     # Description:
##     #   Get data slot from a 'timeSeries' object
##     
##     # Arguments:
##     #   x - a 'timeSeries' object
##     
##     # FUNCTION:
##     
##     # Return Value:
##     series(x)
## }


# ------------------------------------------------------------------------------


## getSeries.timeSeries <- 
##     function(x)
## {
##     # Description:
##     #   Get data slot from a 'timeSeries' object
##     
##     # Arguments:
##     #   x - a 'timeSeries' object
##     
##     # FUNCTION:
##     
##     # Return Value:
##     series(x)
## }


# ------------------------------------------------------------------------------


## "setSeries<-" <-
##     function(x, value)
## {
##     # Description:
##     #   Set data slot to a 'timeSeries' object
## 
##     # Arguments:
##     #   x - a 'timeSeries' object
##     
##     # FUNCTION:
##     
##     # Assign Series Slot:
##     series(x) <- value
##     
##     # Return Value:
##     x
## }


################################################################################
# DEPRECATED


seriesData <-
function(object)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #    Returns the series Data of an ordered data object.

    # Arguments:
    #   object - a 'timeSeries' object

    # Value:
    #    Returns an object of class 'matrix'.

    # FUNCTION:

    # Test:
    if (class(object) != "timeSeries")
        stop("Object is not a time Series")

    # Deprecated
    .Deprecated(new = "series", package = "timeSeries")

    # Get Data Slot:
    ans <- as.matrix(object)

    # Return Value:
    ans
}


################################################################################

