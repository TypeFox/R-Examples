
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
#  holidayTSX                Returns holidays for TSX calendar
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
holidayTSX <-
    function (year = getRmetricsOptions("currentYear"))
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   TSX Holiday Calendar

    # In Canada, the first Monday in August is generally a holiday but it
    # is known by different names in different areas. In Rmetrics it is
    # called CACivicProvincialHoliday()

    # TSX Holidays:
    # http://www.tsx.com/en/market_activity/market_hours.html
    #
    #    * 2007:
    #    * New Year's Day - January 1, 2007
    #    * Good Friday - April 6, 2007
    #    * Victoria Day - May 21, 2007
    #    * Canada Day - July 2, 2007 (for July 1 holiday)
    #    * Civic Day - August 6, 2007
    #    * Labour Day - September 3, 2007
    #    * Thanksgiving Day - October 8, 2007
    #    * Christmas Eve - markets close at 1:00 p.m. ET
    #    * Christmas Day - December 25, 2007
    #    * Boxing Day - December 26, 2007
    #
    #    * 2008:
    #    * New Year's Day - January 1, 2008
    #    * Family Day - February 18, 2008
    #    * Good Friday - March 21, 2008
    #    * Victoria Day - May 19, 2008
    #    * Canada Day - July 1, 2008
    #    * Civic Day - August 4, 2008
    #    * Labour Day - September 1, 2008
    #    * Thanksgiving Day - October 13, 2008
    #    * Christmas Day - December 25, 2008
    #    * Boxing Day - December 26, 2008

    # Trading Hours:
    #   Toronto Stock Exchange and TSX Venture Exchange have trading hours
    #   of 9:30 a.m. to 4:00 p.m. ET, Monday to Friday, with the exception
    #   of the stock market holidays listed below. There is also an extended
    #   session for market participants (Participating Organizations and Members)
    #   from 4:15 to 5:00 p.m. ET each trading day.

    # FUNCTION:

    # Holidays - Years before 2007 are not checked out ...
    holidays = c(
        NewYearsDay(year),
        GoodFriday(year),
        CAVictoriaDay(year),
        CACanadaDay(year),
        CACivicProvincialHoliday(year),
        CAThanksgivingDay(year),
        ChristmasDay(year),
        BoxingDay(year))
    for (y in year)
        if (y >= 2008) holidays = c(holidays, CAFamilyDay(year))
    holidays = sort(holidays)

    # Holidays falling on Saturdays and Sundays:
    holidays =  holidays + (1-isWeekday(holidays))*24*3600
    holidays =  holidays + (1-isWeekday(holidays))*24*3600

    # Add Financial Center:
    holidays <- timeDate(format(holidays),
                         zone = "Toronto", FinCenter = "Toronto")

    # Return Value:
    holidays
}


################################################################################

