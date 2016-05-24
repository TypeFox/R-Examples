
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
#  holiday                   Returns a holiday date of G7 and CH
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
holiday <-
    function(year = getRmetricsOptions("currentYear"), Holiday = "Easter")
{
    # A function implemented by Diethelm Wuertz

    # Description:
    #   Returns the date of a holiday, year may be a vector.

    # Arguments:
    #   year - an integer variable or vector for the year(s) ISO-8601
    #       formatted as "CCYY" as integers.
    #   holiday - a character string naming the holiday. By default
    #       "Easter". Allowable names are the holidays in the G7
    #       countries and Switzerland.

    # Value:
    #   Returns the date of a listed holiday for the selected
    #   "year"(s), an object of class 'timeDate'.

    # Example:
    #   holiday()
    #   holiday(2000:2009, "USLaborDay")
    #   class(holiday())

    # List of Valid Holiday Character Strings:
    #   The following ecclestial and public holidays in
    #       the G7 countries and Switzerland are available:
    #   Holidays Related to Easter:
    #       Septuagesima, Quinquagesima, AshWednesday, PalmSunday,
    #       GoodFriday,  EasterSunday, Easter, EasterMonday,
    #       RogationSunday, Ascension, Pentecost, PentecostMonday,
    #       TrinitySunday CorpusChristi.
    #   Holidays Related to Christmas:
    #       ChristTheKing, Advent1st, Advent1st, Advent3rd,
    #       Advent4th, ChristmasEve, ChristmasDay, BoxingDay,
    #       NewYearsDay.
    #   Other Ecclestical Feasts:
    #       SolemnityOfMary, Epiphany, PresentationOfLord,
    #       Annunciation, TransfigurationOfLord, AssumptionOfMary,
    #       AssumptionOfMary, BirthOfVirginMary, CelebrationOfHolyCross,
    #       MassOfArchangels, AllSaints, AllSouls.
    #   CHZurich - Public Holidays:
    #       CHBerchtoldsDay, CHSechselaeuten, CHAscension,
    #       CHConfederationDay, CHKnabenschiessen.
    #   GBLondon - Public Holidays:
    #       GBMayDay, GBBankHoliday, GBSummerBankHoliday,
    #       GBNewYearsEve.
    #   DEFrankfurt - Public Holidays:
    #       DEAscension, DECorpusChristi, DEGermanUnity, DEChristmasEve,
    #       DENewYearsEve.
    #   FRParis - Public Holidays:
    #       FRFetDeLaVictoire1945, FRAscension, FRBastilleDay,
    #       FRAssumptionVirginMary, FRAllSaints, FRArmisticeDay.
    #   ITMilano - Public Holidays:
    #       ITEpiphany, ITLiberationDay, ITRepublicAnniversary,
    #       ITAssumptionOfVirginMary, ITAllSaints, ITWWIVictoryAnniversary,
    #       ITStAmrose, ITImmaculateConception.
    #   USNewYork/USChicago - Public Holidays:
    #       USNewYearsDay, USInaugurationDay, USMLKingsBirthday,
    #       USLincolnsBirthday, USWashingtonsBirthday, USMemorialDay,
    #       USIndependenceDay, USLaborDay,  USColumbusDay, USElectionDay,
    #       USVeteransDay, USThanksgivingDay, USChristmasDay,
    #       USCPulaskisBirthday, USGoodFriday.
    #   CAToronto/CAMontreal - Public Holidays:
    #       CAVictoriaDay, CACanadaDay, CACivicProvincialHoliday,
    #       CALabourDay, CAThanksgivingDay, CaRemembranceDay.
    #   JPTokyo/JPOsaka - Public Holidays:
    #       JPNewYearsDay, JPGantan, JPBankHolidayJan2, JPBankHolidayJan3,
    #       JPComingOfAgeDay, JPSeijinNoHi, JPNatFoundationDay,
    #       JPKenkokuKinenNoHi, JPGreeneryDay, JPMidoriNoHi,
    #       JPConstitutionDay, JPKenpouKinenBi, JPNationHoliday,
    #       JPKokuminNoKyujitu, JPChildrensDay, JPKodomoNoHi,
    #       JPMarineDay, JPUmiNoHi, JPRespectForTheAgedDay,
    #       JPKeirouNoHi, JPAutumnalEquinox, JPShuubun-no-hi,
    #       JPHealthandSportsDay, JPTaiikuNoHi, JPNationalCultureDay,
    #       JPBunkaNoHi, JPThanksgivingDay, JPKinrouKanshaNohi,
    #       JPKinrou-kansha-no-hi, JPEmperorsBirthday,
    #       JPTennou-tanjyou-bi, JPTennou-tanjyou-bi.
    #   All the holiday functions are listed in the data file "holidays.R"
    #   Additional holidays, which are not yet available there, can be added
    #   to this data base file.

    # FUNCTION:

    # Determine Function:
    nHolidays = length(Holiday)
    ans = NULL
    for (i in 1:nHolidays) {
        FUN = match.fun(Holiday[i])
        ans = c(ans, as.character(FUN(year)))
    }

    # Classify as simple integer ISO date format CCYYMMDD
    ans = timeDate(ans)

    # Return Value:
    ans
}


################################################################################

