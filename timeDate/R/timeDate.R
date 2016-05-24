
# This R package is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This R package is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  timeDate              Creates a 'timeDate' object from given dates
#  setMethod("timeDate", "character",
#  setMethod("timeDate", "timeDate",
#  setMethod("timeDate", "POSIXt",
#  setMethod("timeDate", "Date",
#  setMethod("timeDate", "numeric",
#  setMethod("timeDate", "missing",
#  setMethod("timeDate", "ANY",
# FUNCTION:
#  .formatFinCenter      Internal called by timeDate
#  strptimeDate          Creates for character time stamps a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setGeneric("timeDate",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # Description:
    #   Creates a "timeDate' object from a character vector

    # Arguments:
    #   charvec - a character vector of dates and times. Alternatively
    #       it may be a 'timeDate', a 'Date', or a 'POSIXt' object. In
    #       these cases the argument will be coerced into a character
    #       string or character vector.
    #   format - the format specification of the input character
    #       vector. If set to NULL autodetection will be tried.
    #   zone - the time zone or financial center where the data
    #       were recorded.
    #   FinCenter - a character string with the the location of
    #       the financial center named as "continent/city" where the
    #       data will be used.

    # Value:
    #   Returns a S4 object of class 'timeDate'.

    # Note:
    #   Changeover DST not yet fully implemented!

    # Examples:
    #   timeDate("2004-01-01")
    #   timeDate(c("2004-01-01", "2004-01-01"))
    #   timeDate("2004-01-01 00:00:00")
    #   timeDate("20040101")
    #   timeDate("200401011600")
    #   timeDate("20040101000000")
    #   timeDate("1/1/2004") # American format
    #   timeDate("2004-01-01", FinCenter = "GMT")
    #   timeDate("20040101", FinCenter = "GMT")
    #   td = timeDate("2004-01-01", FinCenter = "GMT"); timeDate(td)
    #   td = timeDate("20040101", FinCenter = "GMT"); timeDate(td)


    standardGeneric("timeDate")
}
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("timeDate", "character",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # Settings and Checks:
    if (zone == "")
        zone <- getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter <- getRmetricsOptions("myFinCenter")

    # ISO Date/Time Format:
    isoDate   <- "%Y-%m-%d"
    isoFormat <- "%Y-%m-%d %H:%M:%S"

    # Autodetect Format :
    if (is.null(format))
        format <- whichFormat(charvec[1])
    if (format %in% c("unknown", "counts")) #-> "counts" catch potential problems from timeSeries
        return(timeDate(NA, zone = zone, FinCenter = FinCenter))

    # Midnight Standard & conversion to isoFormat:
    ct <- midnightStandard2(charvec, format)

    ## Do conversion
    ## YC: .formatFinCenterNum faster than .formatFinCenter
    ## TS: using zone is correct (charvec is converted to GMT)
    num <- .formatFinCenterNum(unclass(ct), zone, type = "any2gmt")

    ## Manually create the POSIXct object:
    ## it is important to set manually the tzone flag,
    num <-
        if (getRversion() >= "2.12.0")
            .POSIXct(num, "GMT")
        else
            structure(num, class = c("POSIXt", "POSIXct"), tzone = "GMT")

    new("timeDate",
        Data = num,
        # Note format is automatically created in
        # initialize,timeDate-method
        FinCenter = as.character(FinCenter))
}
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("timeDate", "timeDate",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # Description:
    #   timeDate

    # if zone not provided, change only the FinCenter in charvec (timeDate)
    if (zone == "") {
        if (FinCenter != "")
            finCenter(charvec) <- FinCenter
        charvec
    } else {
        callGeneric(format(charvec), zone = zone, FinCenter = FinCenter)
    }
}
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("timeDate", "POSIXt",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # Description:
    #   POSIXt

    if (!(zone %in% c("", "GMT", "UTC"))) {
        callGeneric(format(charvec), zone = zone, FinCenter = FinCenter)
    } else {

        # Since zone is not provided consider that charvec is in GMT
        charvec <- as.POSIXct(charvec)
        attr(charvec, "tzone") <- "GMT"

        # FinCenter
        if (FinCenter == "")
            FinCenter = getRmetricsOptions("myFinCenter")

        new("timeDate",
            Data = charvec,
            # Note format is automatically created in
            # initialize,timeDate-method
            FinCenter = as.character(FinCenter))
    }
}
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("timeDate", "Date",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # Description:
    #   Date

    if (!(zone %in% c("", "GMT", "UTC"))) {
        callGeneric(format(charvec), zone = zone, FinCenter = FinCenter)
    } else {

        # Since zone is not provided consider that charvec is in GMT
        charvec <- as.POSIXct(charvec)
        attr(charvec, "tzone") <- "GMT"

        # FinCenter
        if (FinCenter == "")
            FinCenter = getRmetricsOptions("myFinCenter")

        new("timeDate",
            Data = charvec,
            # Note format is automatically created in
            # initialize,timeDate-method
            FinCenter = as.character(FinCenter))
    }
})


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("timeDate", "numeric",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # Description:
    #   numeric

    # DW: Modification of setMethod("timeDate", "numeric") to handle
    #   decimal like inputs (exactly that what "yearmon" does)

    if (!is.null(format) & (format == "%Y" || format == "yearmon" )) {
        # DW: Handels what is known as yearmon format
        # Example:     timeDate(2008+seq(0, 23, by = 1)/12, "yearmon")
        #   Quarterly: timeDate(2008+seq(2, 23, by = 3)/12, format = "%Y")
        # The next 4 lines are borrowed from Zeileis' yearmon()
        year <- floor(charvec + 0.001)
        month <- floor(12 * (charvec - year) + 1 + 0.5 + 0.001)
        dd.start <- as.Date(paste(year, month, 1, sep = "-"))
        # here we concentrate to the end of month date ...
        dd.end <- dd.start + 32 - as.numeric(format(dd.start + 32, "%d"))
        charvec <- as.POSIXct(dd.end, origin = "1970-01-01", tz = "GMT")
    } else {
        charvec <- as.POSIXct(charvec, origin = "1970-01-01", tz = "GMT")
    }

    callGeneric()
}
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("timeDate", "missing",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # Description:
    #   missing

    callGeneric(Sys.time(), format, zone, FinCenter)
}
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("timeDate", "ANY",
    function(charvec, format = NULL, zone = "", FinCenter = "")
{
    # Description:
    #   ANY

    callGeneric(as.character(charvec), format, zone, FinCenter)
}
)


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
# ---------------------------------------------------------------------------- #
.formatFinCenterNum <-
function(num, FinCenter, type = c("gmt2any", "any2gmt"))
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Internal function used by function timeDate()

    if (FinCenter == "GMT" || FinCenter == "UTC")
        return(num)

    type <- match.arg(type)
    signum <- switch(type,
                     "gmt2any" = +1,
                     "any2gmt" = -1)
    ##  otherwise give error

    # Get the DST list from the database:
    try <- try(dst.list <- rulesFinCenter(FinCenter), silent = TRUE)
    if (inherits(try, "try-error"))
        stop(gettextf("'%s' is not a valid FinCenter.", FinCenter))

    offSetIdx <- findInterval(num, dst.list$numeric)
    # consider first DST rule if event occured before
    offSetIdx[offSetIdx < 1] <- 1
    num + signum * dst.list$offSet[offSetIdx]
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.formatFinCenter <-
function(charvec, FinCenter, type = c("gmt2any", "any2gmt"))
{
    # A function implemented by Diethelm Wuertz
    #   thanks to contributions from Martin Maechler

    # Description:
    #   Internal function used by function timeDate()

    if (FinCenter == "GMT" || FinCenter == "UTC")
        return(charvec)

    ## else start working:
    type <- match.arg(type)
    signum <- switch(type,
                     "gmt2any" = +1,
                     "any2gmt" = -1)
    ##  otherwise give error


    # Get the DST list from the database:
    try <- try(dst.list <- rulesFinCenter(FinCenter), silent = TRUE)
    if (inherits(try, "try-error"))
        stop(gettextf("'%s' is not a valid FinCenter.", FinCenter))
    # Update list with last entry:
    z = as.matrix(dst.list)
    z[dim(z)[1], ]
    vec1 = as.vector(c(z[, 1], "2099-01-01 00:00:00"))
    vec2 = as.vector(c(z[, 2], rev(z[, 2])[1]))
    dst.list = data.frame(ruleChanges = as.character(vec1),
    offSet = as.integer(vec2))
    # Extract the dates when DST was changed:
    dst.dates = as.character(dst.list[, 1])
    # Extract the Offsets to GMT
    dst.offsets = as.character(dst.list[, 2])
    # The new dates ar the charvec's:
    new.dates = charvec
    # The new offsets are still unknown:
    new.offsets = rep(NA, length(charvec))
    # Combine all Dates and Offsets:
    dates = c(dst.dates, new.dates)
    offsets = c(dst.offsets, new.offsets)
    # The number of Offsets:
    n = length(dates)
    # Order the Dates:
    o = order(dates)
    # Dates and Offsets in the right order:
    o.dates = dates[o]
    o.offsets = offsets[o]
    # The points at which we have to determine the offsets
    xout = (1:n)[is.na(o.offsets)]
    # The date indexes:
    x = (1:n)[-xout]
    # The corresponding offsets
    y = o.offsets[x]
    # The new offsets:
    yout = approx(x, y , xout, method = "constant")$y
    # All dates:
    m = length(dst.dates)
    # Put them in the right order:
    # Added DW: 2005-05-27
    idx = order(o[which(o>m)])
    offSets = yout[idx]
    dt = strptime(charvec, "%Y-%m-%d %H:%M:%S", tz = "GMT")

    ## Return Value:
    format(dt + signum * offSets, format="%Y-%m-%d %H:%M:%S")
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
strptimeDate <-
function(x, format = whichFormat(x), tz = "")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates for character time stamps a 'timeDate' object

    # Example:
    #   timeDate(); strptimeDate(as.character(Sys.timeDate()))

    # Note:
    #   This function works like strptime.

    # FUNCTION:

    # Check Arguments:
    stopifnot(is.character(x))

    # Settings and Checks:
    if (tz == "")
        tz = getRmetricsOptions("myFinCenter")

    # Create 'timeDate':
    ans = timeDate(x, format, zone = tz, FinCenter = tz)

    # Return Value:
    ans
}


################################################################################

