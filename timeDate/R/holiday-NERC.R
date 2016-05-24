
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
#  holidayNERC               Returns holidays for full-day NERC calendar
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
holidayNERC <-
function (year = getRmetricsOptions("currentYear"), FinCenter = "Eastern")
{
    # A function contributed by Joe W. Byers
    # and modified by Yohan Chalabi

    # Description:
    #   North American Energy Reliabibity Council Holidays calendar
    #   for determining days when the 16 hour on peak hours of a day
    #  for the region are recognized as off peak hours.

    # Details:
    #   Contributed by Joe W Byers:
    #
    #   I have created a calendar function to build out the NERC holidays for
    #   power prices.  NERC: North American Reliability Council.  I started
    #   this because I needed a NYMEX calendar as the text at the bottom
    #   explains, but I am still waiting on NYMEX to reply on my questions.
    #
    #   The NERC holiday function is based on the holidayNYSE function from
    #   Rmetrics and the NERC website rules found at
    #   http://www.nerc.com/~oc/offpeaks.html
    #
    #   I would appreciate any help in verifying this function.  It returns
    #   the correct holidays as shown on their website.  I am not sure if the
    #   historical holidays are correct, especially prior to 1971 when Memorial
    #   day was on 5/30.  I have sent the NERC contact an email with this
    #   question.
    #
    #   I did add a new argument to the function for the FinCenter of the
    #   timedate class because NERC covers all North America.  User will
    #   want to set this to the FinCenter appropriate for their use (Chicago,
    #   Pacific,...).  It would be cool to add NERC regions to the timedate
    #   zone or fincenter lists, but this is another discussion.
    #
    #   The only difference is the holidays
    #   for Independence day and thanks giving.
    #   Independence Day
    #       Monday, July 3*
    #       Tuesday, July 4
    #       (Electronic trading closed Sunday and Monday, July 2 and 3; reopens
    #       6:00 PM, July 4)
    #   Thanksgiving
    #       Thursday, November 23
    #       Friday, November 24
    #       (NYMEX ClearPort(R) and CME Globex(R) open both days)
    #
    #   NYMEX is closed on Monday July 3rd this year since the 4th is on a
    #   Tuesday and Friday the day following Thanksgiving.  This will end in
    #   2007 according to NYMEX.  The Independence day holiday I think will
    #   observe Friday as a holiday when the 4th is on a Thurs but, according
    #   to NYMEX they are not sure for 2007.  Also the webpage
    #   http://www.nymex.com/holida_schedu.aspx notes that Christmas eve will
    #   become a holiday in 2007 but expirations will not change. I am sending
    #   their customer service a request to post a document with the holiday
    #   rules as well as these tables for those of us who worry about a
    #   derivatives calendar.

    # FUNCTION:

    # NERC Holidays:
    holidays <- NULL
    for (y in year) {
        if (y >= 1885)
            holidays <- c(holidays, as.character(USNewYearsDay(y)))
        if (y >= 1885)
            holidays <- c(holidays, as.character(USIndependenceDay(y)))
        if (y >= 1885)
            holidays <- c(holidays, as.character(USThanksgivingDay(y)))
        if (y >= 1885)
            holidays <- c(holidays, as.character(USChristmasDay(y)))
        if (y >= 1887)
            holidays <- c(holidays, as.character(USLaborDay(y)))
        if (y <= 1970)
            holidays <- c(holidays,
        as.character(USDecorationMemorialDay(y)))
        if (y >= 1971)
            holidays <- c(holidays, as.character(USMemorialDay(y)))
    }

    # Sort and Convert to 'timeDate':
    holidays <- sort(holidays)
    ans <- timeDate(holidays, zone = FinCenter, FinCenter = FinCenter)

    # Move Sunday Holidays to Monday:
    posix1 <- as.POSIXlt(ans, tz = "GMT")
    ans <- ans + as.integer(posix1$wday == 0) * 24 * 3600

    # After July 3, 1959, move Saturday holidays to Friday
    # ... except if at the end of monthly/yearly accounting period
    # this is the last business day of a month.
    posix2 <- as.POSIXlt(as.POSIXct(ans, tz = "GMT") - 24 * 3600)
    y <- posix2$year + 1900
    m <- posix2$mon + 1
    calendar <- timeCalendar(y = y + (m + 1) %/% 13,
                             m = m + 1 - (m + 1) %/% 13 * 12,
                             d = 1,
                             zone = "GMT", FinCenter = "GMT")
    lastday <- as.POSIXlt(calendar - 24 * 3600, tz = "GMT")$mday
    lon <- .last.of.nday(year = y, month = m, lastday = lastday, nday = 5)
    ExceptOnLastFriday <- timeDate(format(lon), zone = FinCenter,
                                   FinCenter = FinCenter)
    ans <- ans - as.integer(ans >= timeDate("1959-07-03",
                            zone ="GMT", FinCenter = "GMT") &
                            as.POSIXlt(ans, tz = "GMT")$wday == 0 &
                            ans != ExceptOnLastFriday) * 24 * 3600

    # Remove Remaining Sunday Dates:
    posix3 <- as.POSIXlt(ans, tz = "GMT")
    ans <- ans[!(posix3$wday == 0)]

    # Return Value:
    ans
}


################################################################################

