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
# FUNCTION:              DESCRIPTION:
#  applySeries           Applies a function to blocks of a 'timeSeries'
#  fapply                Applies a function to 'timeSeries' windows
# DEPRECATED:            DESCRIPTION:
#  .applySeries           Applies a function to blocks of a 'timeSeries'
#  .fapply                Applies a function to 'timeSeries' windows
################################################################################


applySeries <-
    function(x, from = NULL, to = NULL, by = c("monthly", "quarterly"),
    FUN = colMeans, units = NULL, format = x@format, zone = x@FinCenter,
    FinCenter = x@FinCenter, recordIDs = data.frame(), title = x@title,
    documentation = x@documentation, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Apply a function to the margins of a 'timeSeries' object

    # Details:
    #   This function can be used to aggregate and coursen a
    #   'timeSeries' object.

    # Arguments:
    #   x - a 'timeSeries' object to be aggregated
    #   from, to - two 'timeDate' position vectors which size the
    #       blocks
    #   by - calendarical block, only active when both 'from'
    #       and 'to' are NULL
    #   FUN - function to be applied, by default 'colMeans'
    #   units - a character vector with column names, allows to
    #       overwrite the column names of the input 'timeSeries'
    #       object.

    # Value:
    #   Returns a S4 object of class 'timeSeries'.

    # Notes:
    #   The size of the 'moving' window and the selection of an
    #   'adj'-acent endpoint are not needed, all the information
    #   is kept in the 'from' and 'to' position vectors.

    # FUNCTION:

    # .Deprecated("aggregate", "timeSeries")

    # Check object:
    if (class(x) != "timeSeries")
        stop("s is not a timeSeries object")

    ###     if (x@format == "counts")
    ###         stop(as.character(match.call())[1],
    ###              " is for time series and not for signal series.")

    # Monthly and Quarterly from and to:
    if (is.null(from) & is.null(to)) {
        if (by[1] == "monthly") {
            # Use monthly blocks:
            from = unique(timeFirstDayInMonth(time(x)))
            to = unique(timeLastDayInMonth(time(x)))
        } else if (by[1] == "quarterly") {
            from = unique(timeFirstDayInQuarter(time(x)))
            to = unique(timeLastDayInQuarter(time(x)))
        } else {
            stop("by must be eiter monthly or quarterly")
        }
        from@FinCenter = to@FinCenter = FinCenter
    }

    # Column Names:
    colNames = units

    # Function:
    fun = match.fun(FUN)

    ###     # Blocks:
    ###     j.pos = as.POSIXct(time(x))
    ###     j.from = as.POSIXct(from)
    ###     j.to = as.POSIXct(to)

    # Blocks:
    j.pos = time(x)
    if (is(j.pos, "timeDate")) {
        j.from = as.timeDate(from)
        j.to = as.timeDate(to)
    } else {
        j.from = as.integer(from)
        j.to = as.integer(to)
    }


    # Iterate:
    pos = time(x)
    rowNames = rownames(x)
    rowBind = NULL
    for (i in seq_len(length(from))) {
        test <- (j.pos >= j.from[i] & j.pos <= j.to[i])
        if (!sum(test)) stop("outsite of range")
        # make sure that cutted is a matrix ...
        cutted = as.matrix(x[test, ])
        # YC : *AND* make sure the matrix is not subbsetted to a vector!!!
        # YC : here it is fine because as.matrix of a timeSeries checks it
        # YC : but prefer to check it one more time at the end of the loop...
        ### if (sum(test)>0) rownames(cutted) <- rowNames[test]
        ans = fun(cutted, ...)
        rowBind = rbind(rowBind, ans)
    }
    stopifnot(NCOL(rowBind) == NCOL(x)) # YC : see above
    # YC : length(to) might not be == NCOL(rowBind)
     if (length(as.character(to)) == NROW(rowBind))
         rownames(rowBind) = as.character(to)

    if (is.null(colNames)) {
        units = x@units
    } else {
        units = colNames }

    # Return Value:
    timeSeries(data = rowBind,  units = units,
        format = format, zone = zone, FinCenter = FinCenter, recordIDs =
        recordIDs, title = title, documentation = documentation, ...)
}


# ------------------------------------------------------------------------------


fapply <-
function(x, from, to, FUN, ...)
{
    # .Deprecated("aggregate", "timeSeries")

    # Check x:
    stopifnot(is(x, "timeSeries"))
    if (x@format == "counts")
        stop(as.character(match.call())[1],
             " is for time series and not for signal series.")

    # Check for missing form/to:
    if(missing(from)) from = start(x)
    if(missing(to)) to = end(x)

    # Return Value:
    applySeries(x = x, from = from, to = to, FUN = FUN, ...)
}


################################################################################
# *** OLD ***
# Check if it is still used somewhere ...


.applySeries <-
    function (x, from = NULL, to = NULL, by = c("monthly", "quarterly"),
    FUN = colMeans, units = NULL, ...)
{
    # Old/Alternative Version

    # Chreck for 'timeSeries' Object:
    stopifnot(is.timeSeries(x),
              is(from, "timeDate") || is.null(from),
              is(to,   "timeDate") || is.null(to))

    # Allow for colMeans:
    if (substitute(FUN) == "colMeans") FUN = "colAvgs"

    # Monthly and Quarterly from and to:
    if (is.null(from) & is.null(to)) {
        by = match.arg(by)
        if (by == "monthly") {
            from = unique(timeFirstDayInMonth(time(x)))
            to = unique(timeLastDayInMonth(time(x)))
        }
        else if (by == "quarterly") {
            from = unique(timeFirstDayInQuarter(time(x)))
            to = unique(timeLastDayInQuarter(time(x)))
        }
        from@FinCenter = to@FinCenter = x@FinCenter
    }

    # Start Cutting Process:
    fun = match.fun(FUN)
    cutted = NULL
    i = 1

    # Find First Interval which is not empty:
    while (is.null(cutted)) {
        cutted = cut(x, from[i], to[i])
        if (!is.null(cutted)) {
            # Non empty Interval:
            ans = fun(cutted, ...)
        }
        i = i + 1
    }
    # Continue up to the end:
    for (j in seq_len(length(from))) {
        cutted = cut(x, from[j], to[j])
        if (!is.null(cutted)) {
            # Non empty Interval:
            newAns = fun(cutted, ...)
            ans = rbind(ans, newAns)
        }
    }

    # Return Value:
    ans
}


################################################################################
# *** OLD ***
# Check if it is still used somewhere ...


.fapply <-
function(x, from, to, FUN, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Applies a function to 'timeSeries' windows

    # Details:
    #   This function can be used to aggregate and coursen a
    #   'timeSeries' object.

    # Arguments:
    #   x - a 'timeSeries' object to be aggregated
    #   from, to - two 'timeDate' position vectors which size the blocks
    #   FUN - function to be applied, by default 'colMeans'

    # Value:
    #   Returns a S4 object of class 'timeSeries' if FUN returns
    #   a time series object, otherwise a list, where the entries
    #   for each window is the output of the function FUN.

    # Notes:
    #   The size of the 'moving' window and the selection of an
    #   'adj'-acent endpoint are not needed, all the information
    #   is kept in the 'from' and 'to' position vectors.

    # FUNCTION:

    # Check object:
    if (class(x) != "timeSeries") stop("s is not a timeSeries object")

    # Monthly and Quarterly from and to:
    if (is.null(from) & is.null(to)) {
        if (by[1] == "monthly") {
            # Use monthly blocks:
            from = unique(timeFirstDayInMonth(time(x)))
            to = unique(timeLastDayInMonth(time(x)))
        } else if (by[1] == "quarterly") {
            from = unique(timeFirstDayInQuarter(time(x)))
            to = unique(timeLastDayInQuarter(time(x)))
        } else {
            stop("by must be eiter monthly or quarterly")
        }
        from@FinCenter = to@FinCenter = x@FinCenter
    }

    # Column Names:
    colNames = units

    # Function:
    fun = match.fun(FUN)

    # Blocks:
    j.pos = as.POSIXct(time(x))
    j.from = as.POSIXct(from)
    j.to = as.POSIXct(to)

    # Iterate:
    y = series(x)
    pos = time(x)
    rowNames = rownames(x)

    # Compute for the first window ...
    i = 1
    test = (j.pos >= j.from[i] & j.pos <= j.to[i])
    # make sure that cutted is a matrix ...
    cutted = as.matrix(y[test, ])
    ### if (sum(test)>0) rownames(cutted) <- rowNames[test]
    ans = fun(cutted, ...)

    if (is.timeSeries(ans)) {
        ## DW can this happen - check ?
        rowBind = ans
        for (i in 2L:length(from)) {
            test = (j.pos >= j.from[1] & j.pos <= j.to[1])
            # make sure that cutted is a matrix ...
            cutted = as.matrix(y[test, ])
            ### if (sum(test)>0) rownames(cutted) <- rowNames[test]
            ans = fun(cutted, ...)
            rowBind = rbind(rowBind, ans)
        }
        rownames(rowBind) = as.character(to)
        if (is.null(colNames)) {
            units = x@units
        } else {
            units = colNames
        }
        # Return Value:
        ans = timeSeries(data = rowBind, charvec = as.character(to),
            units = units, format = format, zone = x@zone, FinCenter =
            x@FinCenter, recordIDs = x@recordIDs, title = x@title,
            documentation = x@documentation, ...)
        return(ans)
    } else {
        listBind = list()
        ## DW [] -> [[]]
        listBind[[1]] = ans
        for (i in 2L:length(from)) {
            test = (j.pos >= j.from[i] & j.pos <= j.to[i])
            # make sure that cutted is a matrix ...
            cutted = as.matrix(y[test, ])
            ### if (sum(test)>0) rownames(cutted) <- rowNames[test]
            ans = fun(cutted, ...)
            ## DW [] -> [[]]
            listBind[[i]] = ans
        }
        # Return Value:
        ans = listBind
        attr(ans, "control") <- list(x = x, from = from, to = to)
        return(invisible(ans))
    }

    # Return Value:
    return()
}


################################################################################
