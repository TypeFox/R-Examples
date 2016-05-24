
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
# You should have received A copy of the GNU Library General
# Public License along with this R package; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                 DATE FUNCTIONS:
#  .fjulian                  Transforms formatted dates to julian day numbers
#  .julian                   Implements SPlus like 'julian'
#  .isPOSIX                  Checks for an object of class POSIX
#  .by2seconds               Converts 'by' string into numeric value of seconds
################################################################################


## .fjulian <-
##     function(fdates, origin = 19600101, order = 'mdy', cc = NULL, swap = 20)
## {
##     # A function implemented by Diethelm Wuertz

##     # Description:
##     #   Transforms formatted dates (fdates) from several formats
##     #   as 8/11/73 11Aug1973, ... into ISO-8601 Gregorian dates
##     #   ... makes use of C-Program char_date.c implemented by
##     #   Terry Therneau

##     # Notes:
##     #   cc - Century, becoming obsolete with the introduction of swap.

##     # Example:
##     #   require(date)
##     #   fdates = c("8/11/73", "08-11-73", "August 11 1973", "Aug11/73")
##     #   .fjulian(fdates)
##     #   fdates = c("11/8/73", "11-08-73", "11 August 1973", "11Aug73")
##     #   .fjulian(fdates, order = 'dmy')

##     # Note:
##     #   Requires R-package "date"

##     # FUNCTION:

##     stopifnot(require("date"))


##     # Formats:
##     order.vec <-
##         switch(order,
##                'ymd'= c(1,2,3),
##                'ydm'= c(1,3,2),
##                'mdy'= c(2,3,1),
##                'myd'= c(2,1,3),
##                'dym'= c(3,1,2),
##                'dmy'= c(3,2,1),
##                stop("Invalid value for 'order' option"))
##     nn = length(fdates)
##     cd <- .C("char_date",
##              as.integer(nn),
##              as.integer(order.vec),
##              as.character(fdates),
##              month = integer(nn),
##              day = integer(nn),
##              year = integer(nn), PACKAGE = "date")[c("month", "day", "year")]

##     yy <- cd$year %% 100

##     # Swap:
##     cc = 19 + trunc(sign(swap-yy)+1)/2
##     year = cc*100 + yy

##     # Origin:
##     cc0 = origin %/% 1000000
##     yymmdd0 = origin - cc0*1000000
##     yy0 = yymmdd0 %/% 10000
##     mm0 = yymmdd0 %/% 100 - yy0*100
##     dd0 = yymmdd0 - yy0*10000 - mm0*100

##     # Result:
##     .julian(cd$month, cd$day, year, origin = c(mm0, dd0, cc0*100+yy0))
## }


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.julian <- function(m, d, y, origin = c(month = 1, day = 1, year = 1960))
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   This function is a synonyme for Splus' "julian()" with the
    #   same list of arguments.

    # Note:
    #   SPlus like function.

    # FUNCTION:

    # Implementation:
    only.origin = all(missing(m), missing(d), missing(y))
    if(only.origin) m = d = y = NULL    # return days since origin
    nms = names(d)
    max.len = max(length(m), length(d), length(y))
    # prepend new origin value and rep out to common max. length:
    m = c(origin[1], rep(m, length = max.len))
    d = c(origin[2], rep(d, length = max.len))
    y = c(origin[3], rep(y, length = max.len))
    # code from julian date in the S book (p.269)
    y = y + ifelse(m > 2, 0, -1)
    m = m + ifelse(m > 2, -3, 9)
    c = y %/% 100
    ya = y - 100 * c
    out = (146097 * c) %/% 4 + (1461 * ya) %/% 4 +
        (153 * m + 2) %/% 5 + d + 1721119
    # now subtract the new origin from all dates
    if(!only.origin) {
        out <- if(all(origin == 0)) out[-1] else out[-1] - out[1]
    }
    names(out) = nms

    # Return Value:
    out
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.isPOSIX <-
    function(x)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Checks for an object of class POSIX

    # FUNCTION:

    # Check:
    ans = inherits(x, "POSIXt")

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.by2seconds <-
    function(by = "1 h")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Convert 'by' string into numeric value of seconds

    # FUNCTION:

    # Convert:
    by = strsplit(by, " ")[[1]]
    byTime = as.integer(by[1])
    byUnits = substr(by[2], 1, 1)
    timeUnits = c(1, 60, 3600)
    names(timeUnits) = c("s", "m", "h")
    bySeconds = byTime * timeUnits[byUnits]
    names(bySeconds) = "secs"

    # Return Value:
    bySeconds
}


################################################################################
