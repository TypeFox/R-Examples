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
# FUNCTION:                     DESCRIPTION:
#  merge,timeSeries,ANY          Merges 'timeSeries' object and ANY
#  merge,timeSeries,missing      Merges 'timeSeries' object and missing
#  merge,timeSeries,timeSeries   Merges two 'timeSeries' objects
#  merge,ANY,timeSeries          Merges ANY and 'timeSeries' object
################################################################################


setMethod("merge", c("timeSeries", "ANY"),
function(x, y, ...)
{
    callGeneric(x, as(y, "timeSeries"), ...)
}
)


# ------------------------------------------------------------------------------


setMethod("merge", c("timeSeries", "missing"),
function(x, y, ...)
{
    x
}
)


# ------------------------------------------------------------------------------


setMethod("merge", c("timeSeries", "numeric"),
function(x, y, ...)
{
    # Deal with names of numeric vectors
    units <-  names(y)
    if (is.null(units))
        units <- paste((substitute(x)), collapse = ".")

    if (length(y) == 1) {
        y <- rep(y, times = nrow(x))
        return(merge(x, timeSeries(y, time(x), units = units), ...))
    } else if (length(y) == nrow(x)) {
        return(merge(x, timeSeries(y, time(x), units = units), ...))
    } else {
        stop("number of rows must match")
    }
}
)


# ------------------------------------------------------------------------------


setMethod("merge", c("timeSeries", "matrix"),
function(x, y, ...)
{
    # deal with names of matrix
    units <- colnames(y)
    if (is.null(units)) {
        units <- paste((substitute(y)), collapse = ".")
        if ((nc <- ncol(y)) > 1)
            units <- paste(units, seq(nc), sep = ".")
    }

    if (nrow(y) != nrow(x))
        stop("number of rows must match")
    else
        merge(x, timeSeries(y, time(x), units = units), ...)
})


# ------------------------------------------------------------------------------'


.merge.timeSeries <- function(x, y, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Merges two objects of class 'timeSeries'

    # Arguments:
    #   x, y - two objects of class 'timeSeries'

    # FUNCTION:
  
    # Compose Attributes - Documentation :
    xAttributes <- getAttributes(x)
    yAttributes <- getAttributes(y)
    Attributes <- .appendList(xAttributes, yAttributes) 
    Documentation <- as.character(date())
    attr(Documentation, "Attributes") <- Attributes

    # Merge:
    if (is.signalSeries(x) | is.signalSeries(y)) {
        data <- merge(getDataPart(x), getDataPart(x))
        return(timeSeries(data = data, units = colnames(data)))
    }

    # Convert to Data Frame
    tx <- as.numeric(time(x), "sec")
    ty <- as.numeric(time(y), "sec")
    df.x <-
        if (prod(dim(rec.x <- x@recordIDs)))
            data.frame(positions = tx, getDataPart(x), rec.x)
        else
            data.frame(positions = tx, getDataPart(x))
    df.y <-
        if (prod(dim(rec.y <- y@recordIDs)))
            data.frame(positions = ty, getDataPart(y), rec.y)
        else
            data.frame(positions = ty, getDataPart(y))

    # Merge as Data Frame:
    df <- merge(df.x, df.y, all = TRUE)
    #-> To avoid problems when using invalid data.frame colnames
    nx <- make.names(colnames(x))
    nxrec <- colnames(rec.x)
    ny <- make.names(colnames(y))
    nyrec <- colnames(rec.y)

    dataIdx <- colnames(df) %in% c(nx, ny)
    recIdx <- colnames(df) %in% c(nxrec, nyrec)

    data <- as.matrix(df[,dataIdx, drop=FALSE])
    recordIDs <- if (any(recIdx)) df[,recIdx, drop=FALSE] else data.frame()
    units <- names(df)[dataIdx]
    charvec <- as.numeric(df[,1])

    # Return Value:
    ans <- timeSeries(data = data, charvec = charvec, units = units,
               zone = "GMT", FinCenter = finCenter(x), recordIDs = recordIDs)
  
    ans@documentation <- Documentation
    ans
}


setMethod("merge", c("timeSeries", "timeSeries"),
          function(x, y, ...) .merge.timeSeries(x, y, ...))

# until UseMethod dispatches S4 methods in 'base' functions
merge.timeSeries <- function(x, y, ...) .merge.timeSeries(x, y, ...)


# ------------------------------------------------------------------------------


setMethod("merge", c("numeric", "timeSeries"),
function(x, y, ...)
{
    # Deal with names of numeric vectors
    units <-  names(x)
    if (is.null(units))
        units <- paste((substitute(x)), collapse = ".")

    if (length(x) == 1) {
        x = rep(x, times = nrow(y))
        return(merge(timeSeries(x, time(y), units = units), y, ...))
    } else if (length(x) == nrow(y)) {
        return(merge(timeSeries(x, time(y), units = units), y, ...))
    } else {
        stop("number of rows must match")
    }
}
)


# ------------------------------------------------------------------------------


setMethod("merge", c("matrix", "timeSeries"),
function(x, y, ...)
{
    # Deal with names of matrix
    units <- colnames(x)
    if (is.null(units)) {
        units <- paste((substitute(x)), collapse = ".")
        if ((nc <- ncol(x)) > 1)
            units <- paste(units, seq(nc), sep = ".")
    }

    if (nrow(x) != nrow(y))
        stop("number of rows must match")
    else
        merge(timeSeries(x, time(y), units = units), y, ...)

})


setMethod("merge", c("ANY", "timeSeries"),
function(x, y, ...)
{
    callGeneric(as(x, "timeSeries"), y, ...)
}
)


################################################################################

