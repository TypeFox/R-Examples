
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
#  holidayZURICH             Returns holidays for ZURICH calendar
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
holidayZURICH <- function(year = getRmetricsOptions("currentYear"))
{
    # A function implemented by Diethelm Wuertz
    # and modified by Yohan Chalabi

    # Description:
    #   Returns a holiday Calendar for Zurich in Switzerland

    # Details:
    #   Inspect the holiday database in "data/holiday.db.R"
    #   ... You can add there additional holidays!
    #       NewYearsDay         Jan, 1st
    #       GoodFriday          2 days before Easter
    #       EasterMonday        1 day after Easter
    #       LaborDay            May, 1st
    #       PentecostMonday     50 days after Easter
    #       ChristmasDay        Dec, 25
    #       BoxingDay           Dec, 26
    #       CHBerchtoldsDay     Jan, 2nd
    #       CHSechselaeuten     3rd Monday in April
    #                           1 week later if it coincides with Easter Monday
    #       CHAscension         39 days after Easter
    #       CHConfederationDay  Aug, 1st
    #       CHKnabenschiessen   2nd Saturday to Monday in Sep

    # Example:
    #   holidayZurich(2008)
    #   holidayZurich(2008)
    #   holidayZurich(2006:2008)

    # FUNCTION:

    # Iterate Years:
    holidays <- c(NewYearsDay(year),
                  GoodFriday(year),
                  EasterMonday(year),
                  LaborDay(year),
                  PentecostMonday(year),
                  ChristmasDay(year),
                  BoxingDay(year),
                  CHBerchtoldsDay(year),
                  CHSechselaeuten(year),
                  CHAscension(year),
                  CHConfederationDay(year),
                  CHKnabenschiessen(year) )

    # Sort and Remove Weekends:
    holidays <- sort(holidays)
    holidays <- holidays[isWeekday(holidays)]

    # Add Financial Center:
    holidays <- timeDate(format(holidays), zone = "Zurich",
                         FinCenter = "Zurich")

    # Return Value:
    holidays
}


################################################################################

