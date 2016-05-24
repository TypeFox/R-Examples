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
# METHOD:                   CREATE A TIMESERIES FROM OTHER OBJECTS:
#  as.timeSeries             Defines method for a 'timeSeries' object
#  as.timeSeries.default     Returns the input
#  as.timeSeries.ts          Transforms a 'data.frame' into a 'timeSeries'
#  as.timeSeries.data.frame  Transforms a 'data.frame' into a 'timeSeries'
#  as.timeSeries.character   Loads and transformas from a demo file
#  as.timeSeries.zoo         Transforms a 'zoo' object into a 'timeSeries'
# METHOD:                   TRANSFORM A TIMESERIES INTO OTHER OBJECTS:
#  as.vector.timeSeries      Converts a univariate 'timeSeries' to a vector
#  as.matrix.timeSeries      Converts a 'timeSeries' to a 'matrix'
#  as.numeric.timeSeries     Converts a 'timeSeries' to a 'numeric'
#  as.data.frame.timeSeries  Converts a 'timeSeries' to a 'data.frame'
#  as.ts.timeSeries          Converts a 'timeSeries' to a 'ts'
#  as.ts.logical             Converts a 'timeSeries' to 'logical'
#  as.list.timeSeries        Converts a 'timeSeries' to 'list'
################################################################################


# YC:
# here keep S3 methods because it should expect an oldClass object as argument


# ------------------------------------------------------------------------------


as.timeSeries <-
    function(x, ...)
{
    UseMethod("as.timeSeries")
}


# ------------------------------------------------------------------------------


as.timeSeries.default <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # FUNCTION:

    # Return Value:
    ans <- timeSeries(x, ...)

    ans
}


setAs("ANY", "timeSeries", function(from) as.timeSeries(from))


# ------------------------------------------------------------------------------


as.timeSeries.ts <-
    function(x, ...)
{
    asTime <- unclass(time(x))
    yearPart <- trunc(asTime)
    decimalPart <- asTime - yearPart
    leapYears <- yearPart%%4 == 0 & (yearPart%%100 != 0 | yearPart%%400 == 0)
    days <- trunc(decimalPart * (365 + leapYears)) + 1

    freq <- frequency(x)
    charvec <-
        if (freq == 4) {
            # Quarterly Data:
            days <- days + 1
            ans <- timeDate(format(strptime(paste(yearPart, days),
                                            format = "%Y %j")),
                            zone = "GMT", FinCenter = "GMT")
            timeLastDayInQuarter(ans)
        } else if (freq == 12) {
            # Monthly Data:
            days <- days + 1
            ans <- timeDate(format(strptime(paste(yearPart, days),
                                            format = "%Y %j")),
                            zone = "GMT", FinCenter = "GMT")
            timeLastDayInMonth(ans)
        } else {
            NA
        }

    # Result:
    tS = timeSeries(x, charvec, ...)
    attr(tS, "ts") <- c(start = round(start(x)),
        frequency = round(frequency(x)), deltat = deltat(x))

    # Return Value:
    tS
}


setAs("ts", "timeSeries", function(from) as.timeSeries(from))


# ------------------------------------------------------------------------------


as.timeSeries.data.frame <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Converts a data.frame into a timeSeries object

    # Notes:
    #   The first column must contain the dates.

    # Examples:
    #   data(bmwRet); head(as.timeSeries(data(bmwRet)))

    # FUNCTION:

    if (all(!(num <- unlist(lapply(x, is.numeric)))))
        stop("x contains no numeric columns")

    # Check if rownames(x) or the first column has a valid ISO-format:
    if (num[1])
        # is.numeric() is better than format == "unkown"
        # which can give wrong result. i.e. whichFormat(0.1253328600)
        suppressWarnings(charvec <- timeDate(rownames(x)))
    else
        suppressWarnings(charvec <- timeDate(as.vector(x[,1])))

    data <- as.matrix(x[, num])
    units <- names(x)[num]
    if (any(!(cl <- num[-1]))) {
        recordIDs <- as.data.frame(x[, !c(TRUE, cl)]) # do not take first column
        names(recordIDs) <- names(x)[!c(TRUE, cl)]
    } else {
        recordIDs <- data.frame()
    }

    # Create Time Series Object:
    timeSeries(data = data,
        charvec = charvec,
        units = units,
        recordIDs = recordIDs, ...)
}


setAs("data.frame", "timeSeries", function(from) as.timeSeries(from))


# ------------------------------------------------------------------------------


as.timeSeries.character <- 
    function(x, ...)
{   # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Example:
    #   as.timeSeries(data(nyse))

    # FUNCTION:

    # Load Demo File - Returns a data frame:
    x <- eval(parse(text = eval(x)))

    # timeSeries:
    ans <- as.timeSeries(x, ...)

    # Return Value:
    ans
}


setAs("character", "timeSeries", function(from) as.timeSeries(from))


# ------------------------------------------------------------------------------


as.timeSeries.zoo <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # FUNCTION:

    # as. timeSeries:

    ans <- timeSeries(data = as.matrix(x),
        charvec = as.character(attr(x, "index")), ...)

    # Return Value:
    ans

}


################################################################################


# YC:
# Since 2.9.0 must define proper S4 methods


.as.matrix.timeSeries <- 
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts a multivariate "timeSeries" to a matrix

    # Arguments:
    #   x - a 'timeSeries' object

    # Value:
    #   Returns the data slot of a 'timesSeries' object as a vector.

    # FUNCTION:

    # Check:
    if (!inherits(x, "timeSeries"))
        stop("x is not a timeSeries object!")

    # Convert:
    ans <- getDataPart(x) # is matrix
    dimnames(ans) <- dimnames(x)

    # Results
    ans
}


setMethod("as.matrix", "timeSeries",
          function(x, ...) .as.matrix.timeSeries(x, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
as.matrix.timeSeries <- function(x, ...) .as.matrix.timeSeries(x, ...)


setAs("timeSeries", "matrix", function(from) as.matrix(from))


# ------------------------------------------------------------------------------


.as.data.frame.timeSeries <- 
    function(x, row.names = NULL, optional = FALSE, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts a multivariate "timeSeries" to a data.frame

    # Arguments:
    #   x - a 'timeSeries' object
    #   row.names, optional - not used

    # Value:
    #   Returns the data slot of a 'timesSeries' object as a data frame.

    # FUNCTION:

    # get rownames from timeSeries
    if (is.null(row.names))
        row.names <- rownames(x)

    if (any(duplicated(row.names)))
        stop("cannot convert to data.frame with duplicate timestamps")

    ans <-
        if (!length(x@recordIDs))
            data.frame(as.list(x), row.names = row.names, ...)
        else
            data.frame(as.list(x), x@recordIDs, row.names = row.names, ...)

    # Return Value:
    ans
}


setMethod("as.data.frame", "timeSeries",
    function(x, row.names = NULL, optional = FALSE, ...)
          .as.data.frame.timeSeries(x, row.names = row.names, optional = optional, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
as.data.frame.timeSeries <- function(x, ...) .as.data.frame.timeSeries(x, ...)


setAs("timeSeries", "data.frame", function(from) as.data.frame(from))


# ------------------------------------------------------------------------------


.as.ts.timeSeries <- 
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Converts a colum from a 'timeSeries' object into an object
    #   of class 'ts'.

    # Example:
    #
    #   x = dummySeries(); as.ts(x)
    #
    #   x = timeSeries(seq(12), timeSequence(by = "month", length.out = 12))
    #   as.ts(x)
    #
    #   x = dummySeries()[c(3,6,9,12),]; as.ts(x)
    #   x = dummySeries()[c(2,5,8,11),]; as.ts(x)
    #   x = dummySeries()[c(1,4,7,10),]; as.ts(x)
    #
    #   x = dummySeries()[c(4,7,10,1),]; as.ts(x)


    # Changes:
    #

    # FUNCTION:

    # check if monthly or quarterly data
    td <- time(x)
    m <- c(timeDate::months(td)) #-> c() to remove attributes
    # (m[1] -1) -> shift vector to match first entry in m
    monthly <- seq(from = m[1]-1, length.out=length(m)) %% 12 + 1
    quarterly <- seq(from = m[1]-1, by = 3, length=length(m)) %% 12 + 1

    # get year of first entry
    y1 <- as.numeric(format(td[1], "%Y"))

    # important to use vector/matrix to avoid troubles with ts()
    data <-
        if (isUnivariate(x))
            as.vector(x)
        else
            matrix(x, ncol = ncol(x))

    if (identical(monthly, m)) # Monthly data
        return(ts(data, start = c(y1, m[1]), frequency = 12))

    if (identical(quarterly, m)) # Quarterly data
        return(ts(data, start = c(y1, m[1]%/%4+1), frequency = 4))

    # other frequencies not implemented yet; return default value
    ans <- ts(data, names = colnames(x))
    attr(ans, "positions") <- time(x)
    ans
}


setMethod("as.ts", "timeSeries", function(x, ...) .as.ts.timeSeries(x, ...))


# until UseMethod dispatches S4 methods in 'base' functions
as.ts.timeSeries <- function(x, ...) .as.ts.timeSeries(x, ...)


setAs("timeSeries", "ts", function(from) as.ts(from))


# ------------------------------------------------------------------------------


# YC: 
# Unneeded since timeSeries inherits from the structure class


# as.logical.timeSeries <- function(x, ...) as.logical(series(x), ...)


# ------------------------------------------------------------------------------


# YC: 
# Important for functions like lapply and sapply to work properly


.as.list.timeSeries <- 
    function(x, ...)
{
    data <- getDataPart(x)
    ncols <- NCOL(data)
    value <- vector("list", ncols)
    for (i in seq.int(ncols)) value[[i]] <- as.vector(data[, i])
    names(value) <- colnames(x)
    value
}


setMethod("as.list", "timeSeries",
          function(x, ...) .as.list.timeSeries(x, ...))

          
# until UseMethod dispatches S4 methods in 'base' functions
as.list.timeSeries <- function(x, ...) .as.list.timeSeries(x, ...)


setAs("timeSeries", "list", function(from) as.list(from))


################################################################################

