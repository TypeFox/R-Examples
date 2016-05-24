
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

# Copyrights (C)
# for this R-port:
#   1999 - Diethelm Wuertz, GPL
#   2007 - Rmetrics Foundation, GPL
#   Diethelm Wuertz <wuertz@phys.ethz.ch>
#   www.rmetrics.org
# for the code accessed (or partly included) from other R-ports:
#   see R's copyright and license files
# for the code accessed (or partly included) from contributed R-ports
# and other sources
#   see Rmetrics's copyright file


################################################################################
# METHODS:             DESCRIPTION:
#  isWeekday            Tests if 'timeDate' falls on a weekday or not
#  isWeekend            Tests if 'timeDate' falls on a weekend or not
#  isBizday             Tests if 'timeDate' falls on a business day or not
#  isHoliday            Tests if 'timeDate' falls on a non-business day or not
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
isWeekday <- 
function(x, wday = 1:5)
{
    # A function implemented by Diethelm Wuertz
    # and improved by Yohan Chalabi

    # Description:
    #   Test if a date is a weekday day or not

    # Arguments:
    #   x - an object of class "timeDate"

    # Value:
    #   returns a logical or a vector of logicals

    # Example:
    #   isWeekday(timeDate("2004-07-01"))
    #   isWeekday(Sys.timeDate())

    # FUNCTION:

    # Test for Weekdays:
    days <- as.POSIXlt(x, tz = "GMT")$wday
    ans <- days %in% wday
    names(ans) <- format(x)

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
isWeekend <- 
function(x, wday = 1:5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Tests if a date is a weekend day or not

    # Arguments:
    #   x - an object of class "timeDate"

    # Value:
    #   Returns a logical or a vector of logicals

    # Example:
    #   isWeekend(timeDate("2004-07-01"))
    #   isWeekend(Sys.timeDate())

    # Changes:
    #

    # FUNCTION:

    # Return Value:
    !isWeekday(x, wday = wday)
}


################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
isBizday <-
function(x, holidays = holidayNYSE(), wday = 1:5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Test if a date is a business day or not

    # Arguments:
    #   x - an object of class "timeDate"
    #   holidays - a holiday calendar

    # Value:
    #   Returns a logical or a vector of logicals

    # Example:
    #   x = timeSequence(from = "2005-05-15", to = "2005-07-15")
    #   h = holiday.NYSE(2005)
    #   cbind(as.character(x), is.bizday(x, h))

    # FUNCTION:

    # Test:
    char.x = substr(as.character(x), 1, 10)
    char.h = substr(as.character(holidays), 1, 10)
    Weekday = as.integer(isWeekday(x, wday = wday))
    nonHoliday = as.integer(!(char.x %in% char.h))

    # Business Days:
    bizdays = as.logical(Weekday*nonHoliday)
    names(bizdays) = x@Data

    # Return Value:
    bizdays
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
isHoliday <-
function(x, holidays = holidayNYSE(), wday = 1:5)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Test if a date is a holiday or not

    # Arguments:
    #   x - an object of class "timeDate"
    #   holidays - a holiday calendar

    # Value:
    #   Returns a logical or a vector of logicals

    # FUNCTION:

    # Return Value:
    !isBizday(x, holidays, wday = wday)
}


################################################################################

