
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
#  Easter                    Returns date of easter or related feasts
# DEPRECATED:               DESCRIPTION:
#  .easterSunday             Easter Algorithm
#  .easter                   Returns date of easter or related feasts
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Easter <-
    function(year = getRmetricsOptions("currentYear"), shift = 0)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns dates of easter or related feasts

    # Arguments:
    #   year - an integer variable or vector for the year(s)
    #       ISO-8601 formatted as "CCYY" where easter or
    #       easter related feasts should be computed.
    #   shift - the number of days shifted from the easter
    #       date. Negative integers are allowed.

    # Value:
    #   Returns the date of Easter shifted by 'shift' days,
    #   ".sdate" formatted, an integer of the form CCYYMMDD.

    # Details:
    #   By default the date of Easter is calculated and returned
    #   in ISO format CCYYMMDD as an integer. Changing shift you
    #   can calculate easter related feasts, e.g. "shift=1" returns
    #   the date of Easter Monday, or "shift=-2" returns the date
    #   of Good Friday.

    # Note:
    #   This algorithm holds for any year in the Gregorian Calendar,
    #   which (of course) means years including and after 1583

    # Examples:
    #   currentYear         # prints current year as integer
    #   .easter()            # date of easter this year
    #   .easter(2000:2009))  # easter for the 2k decade
    #   timeDate(.easter())  # Convert to timeDate
    #   class(.easter())     # what class?

    # Notes:
    #   The variable currentYear is set in ".FirstLib"
    #   Calls ".month.day.year" and ".sjulian"

    # Changes:
    #

    # FUNCTION:

    # Shift and Compute Easter:
    a = year%%19
    b = year%/%100
    c = year%%100
    d = b%/%4
    e = b%%4
    f = (b+8)%/%25
    g = (b-f+1)%/%3
    h = (19*a+b-d-g+15)%%30
    i = c%/%4
    k = c%%4
    l = (32+2*e+2*i-h-k)%%7
    m = (a+11*h+22*l)%/%451
    easter.month = (h+l-7*m+114)%/%31
    p = (h+l-7*m+114)%%31
    easter.day = p+1
    easterSunday = year*10000 + easter.month*100 + easter.day

    # Compose:
    mdy = .month.day.year(.sjulian(easterSunday)+shift)
    ans = as.integer(mdy$year*10000 + mdy$month*100 + mdy$day)

    # Classify as simple integer ISO date format CCYYMMDD
    ans = timeDate(as.character(ans))

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.easterSunday <-
    function(year)
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Computes the 'timeDate' of Easter Sunday

    # FUNCTION:

    # This algorithm holds for any year in the Gregorian Calendar,
    # which (of course) means years including and after 1583
    a = year%%19
    b = year%/%100
    c = year%%100
    d = b%/%4
    e = b%%4
    f = (b+8)%/%25
    g = (b-f+1)%/%3
    h = (19*a+b-d-g+15)%%30
    i = c%/%4
    k = c%%4
    l = (32+2*e+2*i-h-k)%%7
    m = (a+11*h+22*l)%/%451
    easter.month = (h+l-7*m+114)%/%31
    p = (h+l-7*m+114)%%31
    easter.day = p+1

    # Return Value:
    year*10000 + easter.month*100 + easter.day
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
.easter <-
function (year = getRmetricsOptions("currentYear"), shift = 0)
{
    mdy = .month.day.year(.sjulian(.easterSunday(year)) + shift)
    ans = as.integer(mdy$year * 10000 + mdy$month * 100 + mdy$day)
    ans = timeDate(as.character(ans))
    ans
}


################################################################################

