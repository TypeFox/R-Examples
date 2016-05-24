
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
# FUNCTION:                 DESCRIPTION:
#  dayOfWeek                 Returns the day of the week to a 'timeDate' object
#  dayOfYear                 Returns the day of the year to a 'timeDate' object
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
dayOfWeek <-
function(x)
{
    # A function implemented by Diethelm Wuertz
    # and modified by Yohan Chalabi

    # Description:
    #   Returns day of week for time date objects

    # Arguments:
    #   x - an object of class "timeDate"

    # Example:
    #   weekDay(Sys.timeDate())
    #   weekDay(timeSequence("2005-05-15", "2005-07-15"))

    # FUNCTION:
    stopifnot(inherits(x, "timeDate"))

    # Get Day of Week:
    wd <- c("Sun", "Mon", "Tue", "Wed", "Thu", "Fri", "Sat")
    n <- as.POSIXlt(x, tz = "GMT")$wday + 1
    wdays <- wd[n]
    names(wdays) <- format(x)

    # Return Value:
    wdays
}


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
dayOfYear <-
function(x)
{
    # A function implemented by Diethelm Wuertz
    # and modified by Yohan Chalabi

    # Description:
    #   Returns day of week for time date objects

    # Arguments:
    #   x - an object of class "timeDate"

    # FUNCTION:
    stopifnot(inherits(x, "timeDate"))

    # Assign:
    x <- as.POSIXlt(x, tz = "GMT")
    yd <- 1:366
    n <- x$yday + 1
    ydays <- yd[n]
    names(ydays) <- format(x)

    # Return Value:
    ydays
}


################################################################################

