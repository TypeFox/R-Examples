
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
#  ...                       Holiday Dates
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Septuagesima =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, -63)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Quinquagesima =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, -49)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
AshWednesday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, -46)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
PalmSunday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, -7)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
GoodFriday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, -2)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
EasterSunday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
EasterMonday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
RogationSunday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 35)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Ascension =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 39)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Pentecost =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 49)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
PentecostMonday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 50)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
TrinitySunday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 56)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CorpusChristi =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 60)
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ChristTheKing =
function(year = getRmetricsOptions("currentYear")) {
    ans = .on.or.after(year, 11, 20, 0)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Advent1st =
function(year = getRmetricsOptions("currentYear")) {
    ans = .on.or.after(year, 11, 27, 0)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Advent2nd =
function(year = getRmetricsOptions("currentYear")) {
    ans = .on.or.after(year, 12,  4, 0)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Advent3rd =
function(year = getRmetricsOptions("currentYear")) {
    ans = .on.or.after(year, 12, 11, 0)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Advent4th =
function(year = getRmetricsOptions("currentYear")) {
    ans = .on.or.after(year, 12, 18, 0)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ChristmasEve =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1224
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ChristmasDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1225
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
BoxingDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1226
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
SolemnityOfMary =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Epiphany =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0106
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
PresentationOfLord =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0202
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
Annunciation =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0325
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
TransfigurationOfLord =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0806
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
AssumptionOfMary =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0815
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
BirthOfVirginMary =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0908
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CelebrationOfHolyCross =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0914
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
MassOfArchangels =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0929
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
AllSaints =
    function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
AllSouls =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1102
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
NewYearsDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
LaborDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0501
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CHBerchtoldsDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0102
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CHSechselaeuten =
function(year = getRmetricsOptions("currentYear")) {
    ans = NULL
    for (y in year) {
        theDate = .nth.of.nday(y, 4, 1, 3)
        if (as.character(theDate) == as.character(Easter(y, +1))) {
            theDate = .nth.of.nday(y, 4, 1, 4)
        }
        ans = c(ans, theDate)
    }
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CHAscension =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 39)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CHConfederationDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0801
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CHKnabenschiessen =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 9, 1, 2)
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
GBMayDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 5, 1, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
# YC: Note GBBankHoliday returns Spring Bank Holiday
GBBankHoliday =
function(year = getRmetricsOptions("currentYear")) {
    ans = .last.of.nday(year, 5, 31, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
GBSummerBankHoliday =
function(year = getRmetricsOptions("currentYear")) {
    ans = .last.of.nday(year, 8, 31, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
GBMilleniumDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = 19991231
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
DEAscension =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 39)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
DECorpusChristi =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 60)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
DEGermanUnity =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1003
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
DEChristmasEve =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1224
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
DENewYearsEve =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1231
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
FRFetDeLaVictoire1945 =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0508
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
FRAscension =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, 39)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
FRBastilleDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0714
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
FRAssumptionVirginMary =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0815
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
FRAllSaints =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
FRArmisticeDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1111
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ITEpiphany =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0106
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ITLiberationDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0425
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ITAssumptionOfVirginMary =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0815
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ITAllSaints =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ITStAmrose =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1207
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
ITImmaculateConception =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1208
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USNewYearsDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USInaugurationDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0120
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USMLKingsBirthday =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 1, 1, 3)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USLincolnsBirthday =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0212
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USWashingtonsBirthday =
function(year = getRmetricsOptions("currentYear")) {
    ans = year * 10000 + 222
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USMemorialDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .last.of.nday(year, 5, 31, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USIndependenceDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0704
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USLaborDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 9, 1, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USColumbusDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 10, 1, 2)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USElectionDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .on.or.after(year, 11, 2, 2)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USVeteransDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1111
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USThanksgivingDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 11, 4, 4)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USChristmasDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1225
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USCPulaskisBirthday =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 3, 1, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USGoodFriday =
function(year = getRmetricsOptions("currentYear")) {
    ans = Easter(year, -2)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USPresidentsDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 2, 1, 3)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
USDecorationMemorialDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0530
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CAVictoriaDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .on.or.before(year, 5, 24, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CAFamilyDay =
function(year = getRmetricsOptions("currentYear"))
{   # Description:
    #   Adds the new Family Day
    # Note:
    #   Check ...
    #   www.sbhlawyers.com/media/ELD%20Oct%2019%202007%20Public%20Holidays%20and%20Family%20Day.pdf
    #   Family Day will fall on the third Monday of
    #       every February, beginning in 2008.
    # Family Day:
    charvec = paste(year, "02", "01", sep = "-")
    ans = timeNthNdayInMonth(charvec, nday = 1, nth = 3)
    # Return Value:
    ans
}

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CACanadaDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0701
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CACivicProvincialHoliday =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 8, 1, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CALabourDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 9, 1, 1)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CAThanksgivingDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = .nth.of.nday(year, 10, 1, 2)
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
CaRemembranceDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1111
    timeDate(as.character(ans)) }


# ------------------------------------------------------------------------------

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPVernalEquinox <-
function(year = getRmetricsOptions("currentYear"))
{
    # Author:
    #   Parlamis Franklin wrote:
    #   It's me again, with Japanese calendar minutiae I'm sure you've all
    #   been dying to brush up on. The fCalendar functions don't include
    #   the Japanese Vernal Equinox holiday. this is perhaps because there
    #   is no easy way to calculate it. at any rate, here's a function I
    #   wrote to fill the gap.

    # Notes:
    #   Origin and End Date data from
    #   http://aa.usno.navy.mil/data/docs/EarthSeasons.html
    #   The function Vernal.Equinox delivers correct values at the
    #   endpoints of the above data. There may be minor variances
    #   (+/- a few minutes) in the intermediate values, because the
    #   function linearly approximates a phenomenon that is apparently
    #   nonlinear in recorded time.

    Equinox.Origin <- timeCalendar(1992, 3, 20, 8, 48, 0, FinCenter = "GMT")
    Data.EndDate <- timeCalendar(2020, 3, 20, 3, 49, 0, FinCenter = "GMT")
    Total.Seconds <- as.numeric(Data.EndDate-Equinox.Origin)*24*60*60
    Mean.Annual.Seconds <- Total.Seconds / (atoms(Data.EndDate)$Y -
        atoms(Equinox.Origin)$Y)
    Vernal.Equinox <- function(year)
    {
        Equinox.Origin +
        unclass((year-atoms(Equinox.Origin)$Y)*Mean.Annual.Seconds)
    }

    # Nota bene: JP Vernal Equinox is celebrated when the equinox
    #   occurs in the Japanese time zone (see, e.g., 2006, where GMT
    #   Vernal Equinox is on 20 March, but Japanese Equinox holiday is
    #   21 March)

    # Return Value:
    trunc(timeDate(as.character(Vernal.Equinox(year)), FinCenter = "Tokyo"))
}

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPNewYearsDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPGantan =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0101
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPBankHolidayJan2 =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0102
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPBankHolidayJan3 =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0103
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPComingOfAgeDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0115
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPSeijinNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0115
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPNatFoundationDay =
function(year = getRmetricsOptions("currentYear")) {
    ans =year*10000 + 0211
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPKenkokuKinenNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0211
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPGreeneryDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0429
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPMidoriNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0429
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPConstitutionDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0503
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPKenpouKinenBi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0503
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPNationHoliday =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0504
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPKokuminNoKyujitu =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0504
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPChildrensDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0505
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPKodomoNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0505
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPMarineDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0720
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPUmiNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0720
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPRespectForTheAgedDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0915
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPKeirouNOhi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0915
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPAutumnalEquinox =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 0924
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPShuubunNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans =year*10000 + 0924
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPHealthandSportsDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1010
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPTaiikuNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1010
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPNationalCultureDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1103
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPBunkaNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1103
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPThanksgivingDay =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1123
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPKinrouKanshaNoHi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1123
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPEmperorsBirthday =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1123
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPTennouTanjyouBi =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1123
    timeDate(as.character(ans)) }

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
JPBankHolidayDec31 =
function(year = getRmetricsOptions("currentYear")) {
    ans = year*10000 + 1231
    timeDate(as.character(ans)) }


################################################################################

