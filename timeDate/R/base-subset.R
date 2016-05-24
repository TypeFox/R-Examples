
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
# MEHODS:                  SUBSETTING TIMEDATE OBJECTS:
#  setMethod                Extracts/replaces subsets from 'timeDate' objects
#                             signature   missing - missing - ANY
#                             signature   numeric - missing - ANY
#                             signature   logical - missing - ANY
#                             signature   character - missing - ANY
#                             signature   ANY - missing ANY
# "[<-.timeDate"            Extracts/replaces subsets from 'timeDate' objects
# FUNCTION:                DESCRIPTION:
# .subsetCode               Defines codes for different types of subsettings
# .subsetByPython           Subsets a 'timeDate' object by python like indexing
# .subsetBySpan             Subsets a 'timeDate' object by span indexing
################################################################################


# Functions implemented by Yohan Chalabi and Diethelm Wuertz

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("[", signature(x="timeDate", i="missing", j="missing", drop="ANY"),
    function(x, i, j, ..., drop = TRUE) x)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("[", signature(x="timeDate", i="numeric", j="missing", drop="ANY"),
    function(x, i, j, ..., drop = TRUE)
    {
        x@Data <- callGeneric(x@Data, i)
        x
    })


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("[", signature(x="timeDate", i="logical", j="missing", drop="ANY"),
    function(x, i, j, ..., drop = TRUE)
    {
        x@Data <- callGeneric(x@Data, i)
        x
    }
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("[", signature(x="timeDate", i="character", j="missing", drop="ANY"),
    function(x, i, j, ..., drop = TRUE)
    {
        if (length(i) > 1) {
            lt <- lapply(i, function(i, x) "["(x, i), x)
            num <- unlist(lapply(lt, function(td) unclass(td@Data)))
            return(timeDate(num, zone = "GMT", FinCenter = x@FinCenter))
        }
        if (.subsetCode(i) == "SPAN") {
            # Subsetting by Span Indexing:
            return(.subsetBySpan(x, i))
        } else {
            # Subsetting by Python Indexing:
            return(.subsetByPython(x, i))
        }
    }
)


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
setMethod("[", signature(x="timeDate", i="ANY", j="missing", drop="ANY"),
          function(x, i, j, ..., drop = TRUE)
          stop("Not Yet implemented"))


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
"[<-.timeDate" <-
    function(x, ..., value)
{
    # A function implemented by Yohan Chalabi

    # Description:
    #   Extracts or replaces subsets from 'timeDate' objects

    # Arguments:
    #   x - a 'timeDate' object

    # Value:
    #   Returns a subset from a 'timeDate' object.

    # FUNCTION:

    FinCenter <- finCenter(x)

    if (!inherits(value, "timeDate"))
        value <- timeDate(value, zone = FinCenter, FinCenter = FinCenter)

    # Subsets:
    z <- as.POSIXlt(x)
    value <- as.POSIXlt(value)
    val <- "[<-"(z, ..., value)
    val <- as.POSIXct(val)

    # Return Value:
    new("timeDate",
        Data = val,
        format = x@format,
        FinCenter = FinCenter)
}

################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.subsetCode <-
function(subset)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Defines codes for different types of subsettings

    # Details:

    # Python Like Indexing:
    #   Subset:             Code:
    #   ISO8601             00000
    #   ::                  00010
    #   ISO8601::ISO8601    00100
    #   ISO8601::           01000
    #   ::ISO8601           10000

    # Indexing by Spans:
    #   subsets = tolower(c(
    #     "last 1 Month(s)",
    #     "last 1 Week(s)",
    #     "last 1 Day(s)",
    #     "last 1 hour(s)",
    #     "last 1 minute(s)",
    #     "last 1 second(s)"))

    # Example:
    #   .subsetCode("2008-03::")
    #   .subsetCode("last 2 Weeks")

    # Code String:
    if (length(grep("last", subset)) > 0 ) {
        code = "SPAN"
    } else {
        code = paste(
            sign(regexpr("^::[[:digit:]]", subset)[1]+1),
            sign(regexpr("[[:digit:]]::$", subset)[1]+1),
            sign(regexpr("[[:digit:]]::[[:digit:]]", subset)[1]+1),
            as.integer(subset == "::"),
            ## KH : "[a-Z]" is invalid in most locales
            length(grep("[[:alpha:]]", subset)), sep = "")
    }

    # Return Value:
    code
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.subsetByPython <-
function(x = timeCalendar(), subset = "::")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Subsets a 'timeDate' object by python like indexing

    # Arguments:
    #   x - a timeDate object
    #   subset - a python like subset string

    # Example:
    #   .subsetByPython(x, subset = "2008")
    #   .subsetByPython(x, subset = "2008-07")
    #   .subsetByPython(x, subset = "::")
    #   .subsetByPython(x, subset = "2008-07::2008-09")
    #   .subsetByPython(x, subset = "2008-07::")
    #   .subsetByPython(x, subset = "::2008-06")

    # FUNCTION:
    stopifnot(length(subset) == 1)

    # Subset Code:
    code = .subsetCode(subset)

    # Full Vector:
    ans = x

    # Date String:
    date = strsplit(subset, "::")[[1]]

    # 1. DATE
    if(code == "00000") {
        # should return NA if no match found
        idx = grep(date, format(x))
        if (!length(idx))
            ans@Data <- as.POSIXct(NA)
        else
            ans <- x[idx]
    }

    # 2. ::
    if(code == "00010") ans = x

    # Internal Functions:
    .completeStart = function(date) {
        substr(paste0(date, "-01-01"), 1, 10) }
    .completeEnd = function(date) {
        if (nchar(date) == 4)
            paste0(date, "-12-31") else
        if (nchar(date) == 7)
            format(timeLastDayInMonth(paste0(date, "-01"))) else
        if (nchar(date) == 10)
            date }

    # 3. DATE::DATE:
    if(code == "00100")
        ans = window(x, .completeStart(date[1]), .completeEnd(date[2]))

    # 4. DATE::
    if(code == "01000")
        ans = window(x, .completeStart(date[1]), end(x))

    # 5. ::DATE
    if(code == "10000")
        ans = window(x, start(x), .completeEnd(date[2]))

    # Return Value
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.subsetBySpan  <-
function(x = timeCalendar(), subset = "last 3 Months")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Subsets a 'timeDate' object by span indexing

    # Arguments:
    #   x - a timeDate object
    #   subset - a span like subset string

    # Note:
    #   ye[ars]
    #   mo[nths]
    #   da[ys]
    #   ho[urs]
    #   mi[nutes]
    #   se[conds]
    #       ... only "last" spans are implemented

    # Example:
    #   .subsetBySpan(timeCalendar(), "last 2 months")
    #   .subsetBySpan(timeCalendar(), "last 62 days")

    # FUNCTION:
    stopifnot(length(subset) == 1)

    # Get Code:
    code = .subsetCode(subset)
    stopifnot(code == "SPAN")

    # Settings:
    duration = as.numeric(strsplit(subset, " ")[[1]][2])
    len = c(ye = 31622400, mo = 2678400, da = 86400, ho = 3600, mi = 60, se = 1)
    unit = tolower(substr(strsplit(subset, " ")[[1]][3], 1, 2))
    offset = len[unit]*duration

    # Return Value:
    window(x, start = end(x) - offset, end(x))
}


################################################################################

