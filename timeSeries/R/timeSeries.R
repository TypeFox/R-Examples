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
# FUNCTION:                 DESCRIPTION:
#  .signalSeries             Creates a signal series from scratch
#  .timeSeries               Creates a time series from scratch
# METHODS:
#  timeSeries,ANY,ANY
#  timeSeries,matrix,missing
#  timeSeries,matrix,timeDate
#  timeSeries,matrix,numeric
#  timeSeries,matrix,ANY
################################################################################


## .signalSeries : generate units, title, documentation if NULL
##                 data must be a matrix


.signalSeries <-
function(data, charvec, units = NULL, format, zone = "",
    FinCenter = "", recordIDs = data.frame(), title = NULL,
    documentation = NULL, ...)
{
    # Description:

    # Arguments:

    # Note:
    #    it is possible that a ts object is considered as a
    #   matrix when timeSeries method as dispatched. Hence this check

    # FUNCTION:

    if (!is(data, "matrix"))
        data <- as(data, "matrix")

    # Add units, title and Documentation:
    if (is.null(units))
        units <- colnames(data)
    if (is.null(units))
        units <- paste("SS.", seq.int(dim(data)[2]), sep = "")
    if (is.null(title)) 
        title = "Signal Series Object"
    if (is.null(documentation)) 
        documentation = as.character(date())

    # remove rownames of data but keep colnames for
    # functions like var, cov ...
    # Note that if it fails, new("timeSeries" should fail to - normal
    try(dimnames(data) <- list(NULL, units), silent = TRUE)

    ###     new("signalSeries",
    ###         .Data = data,
    ###         units = units,
    ###         recordIDs = recordIDs,
    ###         title = title,
    ###         documentation = documentation)

    new("timeSeries",
        .Data = data,
        units = units,
        positions = numeric(0),
        FinCenter = "",
        format = "counts",
        recordIDs = recordIDs,
        title = title,
        documentation = documentation)
}


# ------------------------------------------------------------------------------


## .timeSeries : generate units, title, documentation if NULL
##               data must be a matrix and charvec a timeDate object


.timeSeries <-
  SERIES <-
function(data, charvec, units = NULL, format, zone = "",
    FinCenter = "", recordIDs = data.frame(), title = NULL,
    documentation = NULL, ...)
{
    # Description:
    #   Creates a time series from scratch
  
    # Arguments:

    # Note:
    #   it is possible that a ts object is considered as a
    #   matrix when timeSeries method as dispatched. Hence this check

    # FUNCTION:

    if (!is(data, "matrix"))
        data <- as(data, "matrix")

    stopifnot(is(charvec, "numeric"))

    # Add units, title and Documentation:
    if (is.null(units))
        units <- colnames(data)
    if (is.null(units))
        units <- paste("TS.", seq.int(dim(data)[2]), sep = "")
    if (is.null(title)) 
        title <- "Time Series Object"
    if (is.null(documentation)) 
        documentation <- as.character(date())
    if (missing(format))
        format <- "%Y-%m-%d"
    if (identical("", FinCenter))
        FinCenter <- "GMT"

    # Remove rownames of data but keep colnames for
    # functions like var, cov ...
    # Note that if it fails, new("timeSeries" should fail to - normal
    try(dimnames(data) <- list(NULL, units), silent = TRUE)

    positions <- charvec # as.numeric(charvec, "sec")
    attributes(positions) <- NULL

    new("timeSeries",
        .Data = data,
        positions =  positions,
        units = units,
        format = format, # charvec@format,
        FinCenter = FinCenter, # charvec@FinCenter,
        recordIDs = recordIDs,
        title = title,
        documentation = documentation)
}


# ------------------------------------------------------------------------------
## missing ANY


setMethod("timeSeries", signature(data = "missing", charvec = "ANY"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
        .signalSeries(data = matrix(NA),
            units = units,
            recordIDs = recordIDs,
            title = title,
            documentation = documentation,
            ...)
    })


# ------------------------------------------------------------------------------
## missing missing


setMethod("timeSeries", signature(data = "missing", charvec = "missing"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
        .signalSeries(data = matrix(NA),
            units = units,
            recordIDs = recordIDs,
            title = title,
            documentation = documentation,
            ...)
    })


# ------------------------------------------------------------------------------
## ANY ANY


setMethod("timeSeries", signature(data = "ANY", charvec = "ANY"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
        data <- as(data, "matrix")
        if (!is(data, "matrix"))
            stop("Could not coerce 'data' to a matrix")
        callGeneric(data = data, charvec = charvec, units = units,
                    format = format, zone = zone, FinCenter =
                    FinCenter, recordIDs = recordIDs, title = title,
                    documentation = documentation, ...)
    })


# ------------------------------------------------------------------------------
## ANY missing


setMethod("timeSeries", signature(data = "ANY", charvec = "missing"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
        data <- as(data, "matrix")
        if (!is(data, "matrix"))
            stop("Could not coerce 'data' to a matrix")
        callGeneric(data = data, units = units,
                    format = format, zone = zone, FinCenter =
                    FinCenter, recordIDs = recordIDs, title = title,
                    documentation = documentation, ...)
    })


# ------------------------------------------------------------------------------
##  matrix missing


setMethod("timeSeries", signature(data = "matrix", charvec = "missing"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
        charvec <- rownames(data)
        if (is.null(charvec)) {
            .signalSeries(data = data, units = units,
                recordIDs = recordIDs, title = title,
                documentation = documentation, ...)
        } else {
            callGeneric(data = data,
                charvec = charvec,
                units = units,
                format = format,
                zone = zone,
                FinCenter = FinCenter,
                recordIDs = recordIDs,
                title = title,
                documentation = documentation,
                ...)
        }
    }
)


# ------------------------------------------------------------------------------
##  matrix timeDate


setMethod("timeSeries", signature(data = "matrix", charvec = "timeDate"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
        if (any(is.na(charvec)))
            return(.signalSeries(data = data, units = units,
                recordIDs = recordIDs, title = title,
                documentation = documentation, ...))
        if (any(!c(zone, FinCenter) %in% ""))
            charvec <- timeDate(charvec, format = format,
                                zone = zone, FinCenter = FinCenter)
        .timeSeries(data = data,
            charvec = as.numeric(charvec, "sec"),
            units = units,
            format = charvec@format,
            FinCenter = charvec@FinCenter,
            recordIDs = recordIDs,
            title = title,
            documentation = documentation, ...)
    }
)


# ------------------------------------------------------------------------------
##  matrix numeric


setMethod("timeSeries", signature(data = "matrix", charvec = "numeric"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
        if (any(!c(zone, FinCenter) %in% "")) {
            td <- timeDate(charvec, zone = zone, FinCenter = FinCenter)
            charvec <- as.numeric(td, "sec")
            FinCenter <- finCenter(td)
        }
        .timeSeries(data = data,
            charvec = charvec,
            units = units,
            FinCenter = FinCenter,
            recordIDs = recordIDs,
            title = title,
            documentation = documentation, ...)
    }
)


# ------------------------------------------------------------------------------
##  matrix      ANY


setMethod("timeSeries", signature(data = "matrix", charvec = "ANY"),
    function (data, charvec, units = NULL, format = NULL, zone = "",
        FinCenter = "", recordIDs = data.frame(), title = NULL,
        documentation = NULL, ...)
    {
    # if charvec NULL returns a signal series
    if (is.null(charvec))
        return(.signalSeries(data = data, units = units,
            recordIDs = recordIDs, title = title,
            documentation = documentation, ...))

    # coerce charvec to timeDate
    charvec <- timeDate(charvec = charvec, format = format,
            zone = zone, FinCenter = FinCenter)

    if (any(is.na(charvec)))
        # Note : there is already a warning in timeDate if there are NA's
        .signalSeries(data = data, units = units,
            recordIDs = recordIDs, title = title,
            documentation = documentation, ...)
    else
        .timeSeries(data = data,
            charvec = as.numeric(charvec, "sec"),
            units = units,
            format = charvec@format,
            FinCenter = charvec@FinCenter,
            recordIDs = recordIDs,
            title = title,
            documentation = documentation,
            ...)
    }
)


################################################################################

