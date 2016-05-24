
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
#  listHolidays              Lists Holidays
# DEPRECATED:               DESCRIPTION:
#  .holidayList              Prints all public and ecclestical holidays
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
listHolidays <-
    function(pattern = ".*")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Prints all public and ecclestical holidays

    # FUNCTION:

    # List:
    holidayList = c(
        "Septuagesima",
        "Quinquagesima",
        "AshWednesday",
        "PalmSunday",
        "GoodFriday",
        "Easter",
        "EasterSunday",
        "EasterMonday",
        "RogationSunday",
        "Ascension",
        "Pentecost",
        "PentecostMonday",
        "TrinitySunday",
        "CorpusChristi",
        "ChristTheKing",
        "Advent1st",
        "Advent2nd",
        "Advent3rd",
        "Advent4th",
        "ChristmasEve",
        "ChristmasDay",
        "BoxingDay",
        "NewYearsDay",
        "SolemnityOfMary",
        "Epiphany",
        "PresentationOfLord",
        "Annunciation",
        "TransfigurationOfLord",
        "AssumptionOfMary",
        "BirthOfVirginMary",
        "CelebrationOfHolyCross",
        "MassOfArchangels",
        "AllSaints",
        "AllSouls",
        "LaborDay",
        "CHBerchtoldsDay",
        "CHSechselaeuten",
        "CHAscension",
        "CHConfederationDay",
        "CHKnabenschiessen",
        "GBMayDay",
        "GBBankHoliday",
        "GBSummerBankHoliday",
        "GBMilleniumDay",
        "DEAscension",
        "DECorpusChristi",
        "DEGermanUnity",
        "DEChristmasEve",
        "DENewYearsEve",
        "FRFetDeLaVictoire1945",
        "FRAscension",
        "FRBastilleDay",
        "FRAssumptionVirginMary",
        "FRAllSaints",
        "FRArmisticeDay",
        "ITEpiphany",
        "ITLiberationDay",
        "ITAssumptionOfVirginMary",
        "ITAllSaints",
        "ITStAmrose",
        "ITImmaculateConception",
        "USDecorationMemorialDay",
        "USPresidentsDay",
        "USNewYearsDay",
        "USInaugurationDay",
        "USMLKingsBirthday",
        "USLincolnsBirthday",
        "USWashingtonsBirthday",
        "USMemorialDay",
        "USIndependenceDay",
        "USLaborDay",
        "USColumbusDay",
        "USElectionDay",
        "USVeteransDay",
        "USThanksgivingDay",
        "USChristmasDay",
        "USCPulaskisBirthday",
        "USGoodFriday",
        "CAVictoriaDay",
        "CACanadaDay",
        "CACivicProvincialHoliday",
        "CALabourDay",
        "CAThanksgivingDay",
        "CaRemembranceDay",
        "JPNewYearsDay",
        "JPGantan",
        "JPBankHolidayJan2",
        "JPBankHolidayJan3",
        "JPComingOfAgeDay",
        "JPSeijinNoHi",
        "JPNatFoundationDay",
        "JPKenkokuKinenNoHi",
        "JPGreeneryDay",
        "JPMidoriNoHi",
        "JPConstitutionDay",
        "JPKenpouKinenBi",
        "JPNationHoliday",
        "JPKokuminNoKyujitu",
        "JPChildrensDay",
        "JPKodomoNoHi",
        "JPMarineDay",
        "JPUmiNoHi",
        "JPRespectForTheAgedDay",
        "JPKeirouNOhi",
        "JPAutumnalEquinox",
        "JPShuubunNoHi",
        "JPHealthandSportsDay",
        "JPTaiikuNoHi",
        "JPNationalCultureDay",
        "JPBunkaNoHi",
        "JPThanksgivingDay",
        "JPKinrouKanshaNoHi",
        "JPEmperorsBirthday",
        "JPTennouTanjyouBi",
        "JPBankHolidayDec31")

    # Financial Centers:
    if (pattern == "*") pattern = "\\\\*"

    # Return Value:
    sort(as.character(holidayList[grep(pattern = pattern, x = holidayList)]))
}


# ------------------------------------------------------------------------------


# Keep for compatibility of older Rmetrics versions ...


.holidayList <- listHolidays()


################################################################################

