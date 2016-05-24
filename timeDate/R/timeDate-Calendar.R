
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
# FUNCTION:                 GENERATION OF TIMEDATE OBJECTS:
#  timeCalendar              Creates a 'timeDate' object from calendar atoms
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
timeCalendar <-
    function(y = getRmetricsOptions("currentYear"),
    m = 1:12, d = 1, h = 0, min = 0, s = 0,
    zone = "", FinCenter = "")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Creates a 'timeDate' object from calendar atoms

    # Arguments:
    #   y - calendar years (e.g. 1997), defaults are 1960.
    #   m - calendar months (1-12), defaults are 1.
    #   d - calendar days (1-31), defaults are 1.
    #   h - hours of the days (0-23), defaults are 0.
    #   min - minutes of the days (0-59), defaults are 0.
    #   s - seconds of the days (0-59), defaults are 0.
    #   FinCenter - a character sting with the the location of the
    #       financial center named as "continent/city"

    # Value:
    #   Returns a 'timeDate' object corresponding to the "atomic"
    #   inputs. For the default arguments the first day in each
    #   month of the current year will be returned.

    # Details:
    #   Creates a 'timeDate' object from date as month, day, year and
    #   time of day as hours, and minutes [seconds, milliseconds]

    # Note:
    #   The 'zone' where the data were recorded is fixed to myFincenter!
    #   The argument list has ISO-8601 ordering!
    #   ms - Milliseconds is not supported.

    # Example:
    #   x = timeCalendar(y = 2000, h = rep(16,12))
    #   x = timeCalendar(m = c(3,4,5), d = c(12,15,7), y = c(1998,1997,2004))
    #   x = timeCalendar(h = c(9,14), min = c(15,23))

    # FUNCTION:

    # Settings and Check:
    trace = FALSE
    if (zone == "")
        zone <- getRmetricsOptions("myFinCenter")
    if (FinCenter == "")
        FinCenter <- getRmetricsOptions("myFinCenter")

    if (is.null(h) & is.null(min) & is.null(s)) zone = FinCenter

    # Check Input:
    len = c(length(m), length(d), length(y), length(h), length(min), length(s))
    data.len = max(len)
    if (data.len < 1)
        stop("No arguments defined!")
    if (any((data.len %% len[len > 0]) != 0))
        stop("Arguments have incompatible lengths")

    # Make All Arguments the Same Length:
    if (len[1] == 0) m = 1
    if (len[2] == 0) d = 1
    if (len[3] == 0) y = 1960
    if (len[4] == 0) h = 0
    if (len[5] == 0) min = 0
    if (len[6] == 0) s = 0

    # Presettings:
    # m = rep(m, length = data.len)
    # d = rep(d, length = data.len)
    # y = rep(y, length = data.len)
    # h = rep(h, length = data.len)
    # min = rep(min, length = data.len)
    # s = rep(s, length = data.len)
    # DW 2006-03-13
    if (length(m) < data.len) m = rep(m, length = data.len)
    if (length(d) < data.len) d = rep(d, length = data.len)
    if (length(y) < data.len) y = rep(y, length = data.len)
    if (length(h) < data.len) h = rep(h, length = data.len)
    if (length(min) < data.len) min = rep(min, length = data.len)
    if (length(s) < data.len) s = rep(s, length = data.len)

    # Date-Time Strings:
    # Note Format is always of type  "%Y%m%d%H%M%S"  !
    CCYYMMDD = as.integer(y*10000 + m*100 + d)
    chardate = as.character(CCYYMMDD)
    xhhmmss = as.integer(1000000 + h*10000 + min*100 + s)

    # Date and Date/Time Checks:
    if (mean(xhhmmss) == 1000000) {
        chartime = substr(as.character(xhhmmss), 2, 7)
        charvec = as.vector(chardate)
        format = "%Y%m%d"
    } else {
        chartime = substr(as.character(xhhmmss), 2, 7)
        charvec = paste0(as.vector(chardate), as.vector(chartime))
        format = "%Y%m%d%H%M%S"
    }

    # Return Value:
    ans <- timeDate(charvec = charvec, format = format,
        zone = zone, FinCenter = FinCenter)

    ans
}


################################################################################

