
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
#  timeLastDayInMonth        Computes the last day in a given month and year
#  timeFirstDayInMonth       Computes the first day in a given month and year
#  timeLastDayInQuarter      Computes the last day in a given quarter and year
#  timeFirstDayInQuarter     Computes the first day in a given quarter and year
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
timeLastDayInMonth <-
    function(charvec, format = "%Y-%m-%d", zone = "",
    FinCenter = "")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the last day in a given month and year

    # Arguments:
    #   charvec - a character vector of dates and times.
    #   format - the format specification of the input character vector.
    #   FinCenter - a character string with the the location of the
    #       financial center named as "continent/city".

    # Value:
    #   Returns the last day in a given month and year as a
    #   'timeDate' object.

    # FUNCTION:
    if (zone == "")
        zone = getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter = getRmetricsOptions("myFinCenter")

    # Last day of month:
    last.day = c(31,28,31, 30,31,30, 31,31,30, 31,30,31)
    lt = strptime(charvec, format, tz = "GMT")
    y = 1900 + lt$year
    leap.year = (y%%4 == 0 & (y%%100 != 0 | y%%400 == 0))
    leap.day = as.integer(leap.year)*as.integer(lt$mon == 1)
    lt$mday = last.day[1 + lt$mon] + leap.day

    # Return Value:
    timeDate(format(lt), format = "%Y-%m-%d", zone = zone, FinCenter = FinCenter)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
timeFirstDayInMonth <-
    function(charvec, format = "%Y-%m-%d", zone = "",
    FinCenter = "")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the last day in a given month and year

    # Changes:
    #

    # FUNCTION:
    if (zone == "")
        zone = getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter = getRmetricsOptions("myFinCenter")

    # First Day In Month:
    lt = strptime(charvec, format, tz = "GMT")
    lt$mday = 1

    # Return Value:
    timeDate(format(lt), format = "%Y-%m-%d", zone = zone, FinCenter = FinCenter)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
timeLastDayInQuarter <-
    function(charvec, format = "%Y-%m-%d", zone = "",
    FinCenter = "")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the last day in a given month and year

    # FUNCTION:
    if (zone == "")
        zone = getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter = getRmetricsOptions("myFinCenter")

    # First Day in Month:
    charvec = timeFirstDayInMonth(charvec = charvec, format = format,
        FinCenter = FinCenter)

    # Last Day in Quarter:
    lt = strptime(charvec, format, tz = "GMT")
    last.quarter = rep(c(3,6,9,12), each = 3) - 1
    lt$mon = last.quarter[1 + lt$mon]
    charvec = timeDate(format(lt), format = "%Y-%m-%d", zone = zone,
        FinCenter = FinCenter)

    # Return Value:
    timeLastDayInMonth(charvec = charvec, format = format,
        zone = zone, FinCenter = FinCenter)
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
timeFirstDayInQuarter <-
    function(charvec, format = "%Y-%m-%d", zone = "",
    FinCenter = "")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the last day in a given month and year

    # Changes:
    #

    # FUNCTION:
    if (zone == "")
        zone = getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter = getRmetricsOptions("myFinCenter")

    # First Day in Month:
    charvec = timeFirstDayInMonth(charvec =charvec, format = format,
        FinCenter = FinCenter)

    # First Day in Quarter:
    lt = strptime(charvec, format, tz = "GMT")
    first.quarter = rep(c(1,4,7,10), each = 3) - 1
    lt$mon = first.quarter[1 + lt$mon]

    # Return Value:
    timeDate(format(lt), format = "%Y-%m-%d", zone = zone, FinCenter = FinCenter)
}


################################################################################

