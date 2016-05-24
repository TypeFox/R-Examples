
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
#  holidayLONDON             Returns holidays for British Bank Holidays
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
holidayLONDON <- function (year = getRmetricsOptions("currentYear")) {

    # function implemented by Menon Murali

    holidays <- NULL
    for (y in year) {
        if (y >= 1834 & y <= 1870) {
            # 1 May, 1 November, Good Friday and Christmas are the
            # only England Bank Holidays
            dts <- c(paste0(y, "-05-01"), paste0(y, "-11-01"))
            holidays <- c(holidays, dts, as.character(GoodFriday(y)),
                          as.character(ChristmasDay(y)))
        }
        if (y >= 1871) {
            # Good Friday and Easter Monday
            holidays <- c(holidays, as.character(GoodFriday(y)),
                          as.character(EasterMonday(y)))
            if (y <= 1964) {
                # Whit Monday, which is exactly 50 days after Easter
                holidays <- c(holidays, as.character(Easter(y, 50)))
                # First Monday in August
                lon <- timeDate(.on.or.after(y, 8, 1, 1), zone = "London",
                                FinCenter = "Europe/London")
                holidays <- c(holidays, as.character(lon))
            } else {
                # Last Monday in May replaces Whit Monday
                if (y == 2002) {
                    # Last Monday in May holiday moved to June 3, and
                    # Queen's Jubilee on June 4
                    dts <- c(paste0(y, "-06-03"),
                             paste0(y, "-06-04"))
                    holidays <- c(holidays, dts)
                } else if (y == 2012) {
                    # Last Monday in May holiday moved to June 4, and
                    # Queen's Diamond Jubilee on June 5
                    dts <- c(paste0(y, "-06-04"),
                             paste0(y, "-06-05"))
                    holidays <- c(holidays, dts)
                } else {
                    lon <- timeDate(.last.of.nday(y, 5, 31, 1), zone = "London",
                                    FinCenter = "Europe/London")
                    holidays <- c(holidays, as.character(lon))
                }

                # Last Monday in August replaces first Monday in August
                lon <- timeDate(.last.of.nday(y, 8, 31, 1), zone = "London",
                                FinCenter = "Europe/London")
                holidays <- c(holidays, as.character(lon))
            }

            # Not entirely sure when Mon/Tue began to be given as
            # holiday when Christmas/Boxing Day fell on
            # Saturday/Sun. I'm assuming this was after the 1971
            # Banking and Financial Dealings Act which established the
            # new holiday schedule.
            if (y < 1970) {
                # Christmas and Boxing Day
                holidays <- c(holidays, as.character(ChristmasDay(y)),
                              as.character(BoxingDay(y)))
            } else {
                posix1 <- as.POSIXlt(ChristmasDay(y))
                # If Christmas on Saturday or Sunday, then the following Monday and Tuesday are holidays
                if (posix1$wday == 0) { # Christmas Sunday
                    holidays <- c(holidays,
                                  as.character(ChristmasDay(y) + (1 : 2) * 86400))
                } else if (posix1$wday == 6) { #Christmas Saturday
                    holidays <- c(holidays,
                                  as.character(ChristmasDay(y) + (2 : 3) * 86400))
                } else if (posix1$wday == 5) {# Christmas Friday
                    # the next Monday is a holiday
                    holidays <- c(holidays,
                                  as.character(ChristmasDay(y) + c(0, 3) * 86400))
                } else {
                    holidays <- c(holidays, as.character(ChristmasDay(y)),
                                  as.character(BoxingDay(y)))
                }
            }

            if (y >= 1974) {
                # New Year's Day: if it falls on Sat/Sun, then is
                # moved to following Monday
                posix1 <- as.POSIXlt(NewYearsDay(y))
                if (posix1$wday == 0 | posix1$wday == 6) {
                    lon <- timeDate(.on.or.after(y, 1, 1, 1), zone = "London",
                                    FinCenter = "Europe/London")
                    holidays <- c(holidays, as.character(lon))
                } else {
                    holidays <- c(holidays, as.character(posix1))
                }
            }
            if (y >= 1978) {
                if (y == 1981) {
                    # Royal wedding was a public holiday
                    dts <- paste0(y, "-07-29")
                    holidays <- c(holidays, dts)
                }
                if (y == 2011) {
                    # Royal wedding declared a public holiday
                    dts <- paste0(y, "-04-29")
                    holidays <- c(holidays, dts)
                }
                # First Monday of May became a bank holiday
                if (y == 1995) {
                    # Was moved to May 8 to celebrate VE Day's 50th anniversary
                    dts <- paste0(y, "-05-08")
                    holidays <- c(holidays, dts)
                } else {
                    lon <- timeDate(.on.or.after(y, 5, 1, 1), zone = "London",
                                    FinCenter = "Europe/London")
                    holidays <- c(holidays, as.character(lon))
                }
            }
        }
    }

    holidays <- sort(holidays)
    ans <- timeDate(format(holidays), zone = "London",
                    FinCenter = "Europe/London")
    posix1 <- as.POSIXlt(ans, tz = "GMT")
    ans[!(posix1$wday == 0 | posix1$wday == 6)] # Remove any Saturdays/Sundays
}
