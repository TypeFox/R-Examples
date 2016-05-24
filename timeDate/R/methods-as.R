
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
# METHOD:                   DESCRIPTION:
#  as.timeDate               Implements Use Method
#  as.timeDate.default       Default Method
#  as.timeDate.POSIXt        Returns a 'POSIXt' object as 'timeDate' object
#  as.timeDate.Date          Returns a 'Date' object as 'timeDate' object
#  as.character.timeDate     Returns a 'timeDate' object as 'character' string
#  as.double.timeDate        Returns a 'timeDate' object as 'numeric' object
#  as.data.frame.timeDate    Returns a 'timeDate' object as 'data.frame' object
#  as.list.timeDate          Returns a 'timeDate' object as 'list' object
#  as.POSIXct.timeDate       Returns a 'timeDate' object as 'POSIXct' object
#  as.POSIXlt.timeDate       Returns a 'timeDate' object as 'POSIXlt' object
#  as.Date.timeDate          Returns a 'timeDate' object as 'Date' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.timeDate <-
    function(x, zone = NULL, FinCenter = NULL)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    UseMethod("as.timeDate")
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.timeDate.default <-
    function(x, zone = "", FinCenter = "")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns default object as 'timeDate' object

    # Arguments:
    #   x - a 'timeDate' object

    # Value:
    #   Returns 'x' as a 'timeDate' object.

    # FUNCTION:

    # as timeDate:
    if (zone == "")
        zone <- getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter <- getRmetricsOptions("myFinCenter")

    # Return Value:
    timeDate(as.character(x), zone = zone, FinCenter = FinCenter)
}


setAs("ANY", "timeDate", function(from) as.timeDate.default(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.timeDate.timeDate <-
    function(x, zone = x@FinCenter, FinCenter = "")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns default object as 'timeDate' object

    # Arguments:
    #   x - a 'timeDate' object

    # Value:
    #   Returns 'x' as a 'timeDate' object.

    stopifnot(is(x, "timeDate"))

    if (FinCenter == "")
        FinCenter <- getRmetricsOptions("myFinCenter")
    if (zone != x@FinCenter)
        warning("argument zone is ignored and FinCenter\n of timeDate is used as zone")
    ## FIXME : and now?   'zone' is *NOT* ignored!

    # Return as timeDate:
    timeDate(as.character(x), zone = zone, FinCenter = FinCenter)
}


# setAs("timeDate", "timeDate", function(from) as.timeDate.timeDate(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.timeDate.Date <- function(x, zone = "", FinCenter = "")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns a 'Date' object as 'timeDate' object

    # Arguments:
    #   x - a 'Date' (or 'POSIXt') object

    # Value:
    #   Returns 'x' as  timeDate object.

    # FUNCTION:
    if (zone == "")
        zone <- getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter <- getRmetricsOptions("myFinCenter")

    # Return as timeDate:
    timeDate(x, zone = zone, FinCenter = FinCenter)
}


setAs("Date", "timeDate", function(from) as.timeDate.Date(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.timeDate.POSIXt <- function(x, zone = "", FinCenter = "")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    timeDate(x, zone = zone, FinCenter = FinCenter)
}


setAs("POSIXt", "timeDate", function(from) as.timeDate.POSIXt(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.character.timeDate <-
function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    format(x, ...)
}


setAs("timeDate", "character", function(from) as.character.timeDate(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.double.timeDate <-
    function(x,
    units = c("auto", "secs", "mins", "hours", "days", "weeks"), ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns a 'timeDate' object as 'numeric' vector

    # Arguments:
    #   x - a 'timeDate' object
    #   units - a character string denoting in which units the
    #       elements of the numeric vector are measured

    # Value:
    #   Returns 'x' as a numeric vector.

    # FUNCTION:

    units <- match.arg(units)

    if (units == "secs") {
        ans <- c(unclass(as.POSIXct(x)))
    } else {
        ct = as.POSIXct(x)
        origin = as.POSIXct("1970-01-01", tz = "GMT")
        dt = difftime(ct, origin, tz = "GMT", units = units)
        units <- attr(dt, "units")
        ans = as.double(dt)
    }

    attr(ans, "FinCenter") <- "GMT"
    attr(ans, "units") <- units
    attr(ans, "origin") <-
        switch(units,
               secs  = "1970-01-01 00:00:00 GMT",
               mins  = "1970-01-01 00:00 GMT",
               hours = "1970-01-01 00:00 GMT",
               days  = "1970-01-01 GMT",
               weeks = "1970-01-01 GMT")

    # Return Value:
    ans
}


setAs("timeDate", "numeric", function(from) as.double.timeDate(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.data.frame.timeDate <-
    function(x, ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns a 'timeDate' object as data frame

    # Arguments:
    #   x - a 'timeDate' object

    # Value:
    #   Returns 'x' as a data frame.

    # FUNCTION:

    # Check Class Type:
    stopifnot(inherits(x, "timeDate"))

    # Data Frame:
    ans <- as.data.frame.POSIXlt(x@Data, ...)
    nm <- paste(deparse(substitute(x), width.cutoff = 500), collapse = " ")
    colnames(ans) <- paste0(x@FinCenter, ":", nm)
    attr(ans, "control") <- c(FinCenter = x@FinCenter)

    # Return Value:
    ans
}


setAs("timeDate", "data.frame", function(from) as.data.frame.timeDate(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.list.timeDate <-
    function(x, ...)
{
    # A function implemented by Yohan Chalabi and Diethelm Wuertz

    # Description:
    #   Returns a 'timeDate' object as list
    # important for functions like sapply and lapply

    # Arguments:
    #   x - a 'timeDate' object

    # Value:
    #   Returns 'x' as a data frame.

    # FUNCTION:

    # Check Class Type:
    stopifnot(inherits(x, "timeDate"))

    ans <- vector("list", length(x))

    for (i in seq(length(x)))
        ans[i] <- x[i]

    # Return Value:
    ans
}


setAs("timeDate", "list", function(from) as.list.timeDate(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.POSIXct.timeDate <-
    function(x, tz = "", ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns a 'timeDate' object as POSIXct object

    # Arguments:
    #   x - a 'timeDate' object
    #   tz - a timezone specification to be used for the conversion.
    #       (If tz is used, the method does not consider the
    #        FinCenter of timeDate)

    # Value:
    #   Returns 'x' as an object of class 'POSIXct'.

    # FUNCTION:

    # Check Class Type:
    if (!inherits(x, "timeDate"))
        stop("Wrong class type")

    # POSIXct:
    FinCenter <- finCenter(x)
    if (identical(tz, "")) {
        ans <- getDataPart(x)
        attr(ans, "control") <- c(FinCenter = finCenter(x))
    } else {
        num <- .formatFinCenterNum(as.numeric(getDataPart(x)),
                                   FinCenter, type = "gmt2any")
        ans <- as.POSIXct(num, origin = "1970-01-01", tz = tz)
        attr(ans, "control") <- c(FinCenter = FinCenter)
    }

    # Return Value:
    ans
}

setAs("timeDate", "POSIXct", function(from) as.POSIXct.timeDate(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.POSIXlt.timeDate  <-
    function(x, tz = "", ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns a 'timeDate' object as 'POSIXlt' object

    # Arguments:
    #   x - a 'timeDate' object
    #   tz - a timezone specification to be used for the conversion.
    #       (If tz is used, the method does not consider the
    #        FinCenter of timeDate)

    # Value:
    #   Returns 'x' as an object of class 'POSIXct'.

    # Note:
    #   be careful if S4 method are defined because arguments of generic
    #   function has changed for as.POSIXlt since R-2.6

    # FUNCTION:

    # Set Timezone to GMT:
    ans = as.POSIXlt(as.POSIXct(x = x, tz = tz))

    # Return Value:
    ans
}


setAs("timeDate", "POSIXlt", function(from) as.POSIXlt.timeDate(from))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
as.Date.timeDate <-
    function(x, method = c("trunc", "round", "next"), ...)
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns a 'timeDate' object as 'Date' object

    # Arguments:
    #   x - a 'timeDate' object
    #   method - a character string denoting the method how to
    #       compute the 'Date' object.

    # Value:
    #   Returns 'x' as an object of class 'POSIXct'.

    # FUNCTION:

    # as Date:
    method = match.arg(method)

    # Note:
    #   # Note: one must be careful when converting to Date with tzone.
    #   td <- timeDate("2008-12-11 00:00:01", zone = "Zurich", FinCenter = "Zurich")
    #   ct <- td@Data
    #   attr(ct, "tzone") <- "Europe/Zurich"
    #   # ct and td should be identical
    #   ct; td
    #   # but
    #   as.Date(ct) # trunc on previous day because trunc in GMT
    #   as.Date(td) # trunc in the current FinCenter !

    ans <- switch(method,
                  trunc = as.Date(format(trunc(x)), ...),
                  round = as.Date(format(round(x)), ...),
                  "next" = as.Date(format(trunc(x)), ...) + 1)

    # Add Attribute:
    attr(ans, "control") <- c(method = method, FinCenter = x@FinCenter)

    # Return Value:
    ans
}


setAs("timeDate", "Date", function(from) as.Date.timeDate(from))


################################################################################

