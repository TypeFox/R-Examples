
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
# FUNCTION:               DESCRIPTION:
#  periods                 Returns start and end dates for rolling periods
#  periodicallyRolling     Returns start and end dates for periodically periods
#  monthlyRolling          Returns start and end dates for monthly periods
################################################################################

# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
periods <-
    function (x, period = "12m", by = "1m", offset = "0d")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns start and end dates for a rolling periods

    # Arguments:
    #   x - an object of class timeDate
    #   period - a span string, consisting of a length integer
    #       and a unit value, e.g. "52w" for 52 weeks.
    #   by - a span string, consisting of a length integer
    #       and a unit value, e.g. "4w" for 4 weeks.
    #   offset - a span string, consisting of a length integer
    #       and a unit value, e.g. "0d" for no offset.

    # Details:
    #   Periodically Rolling - Allowed unit values are "m" for
    #       4 weeks, "w" for weeks, "d" for days, "H" for hours, "M"
    #       for minutes, and "S" for seconds.
    #   Monthly Calendar Rolling - The only allowed allowed unit
    #       value is "m" for monthly periods. Express a quarterly
    #       period by "3m", a semester by "6m", a year by "12m" etc.

    # Example:
    #   x = time(as.timeSeries(data(smallcap.ts)))
    #   periods(x, "12m", "1m")
    #   periods(x, "52w", "4w")

    # FUNCTION:

    # Check x:
    stopifnot(is(x, "timeDate"))

    # Check Periods:
    Names = c("m", "w", "d", "H", "M", "S")
    periodUnit = gsub("[ 0-9]", "", period, perl = TRUE)
    stopifnot(periodUnit %in% Names)
    offsetUnit = gsub("[ 0-9]", "", offset, perl = TRUE)
    stopifnot(offsetUnit %in% Names)
    byUnit = gsub("[ 0-9]", "", by, perl = TRUE)
    stopifnot(byUnit %in% Names)

    # Rolling Periods:
    if (periodUnit == "m" & byUnit == "m") {
        ans = monthlyRolling(x, period, by)
    } else {
        ans = periodicallyRolling(x, period, by)
    }

    # Return Value:
    ans
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
periodicallyRolling <-
    function(x, period = "52w", by = "4w", offset = "0d")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns start and end dates for a rolling periods

    # Arguments:
    #   x - an object of class timeDate
    #   period - a span string, consisting of a length integer
    #       and a unit value, e.g. "52w" for 52 weeks.
    #   by - a span string, consisting of a length integer
    #       and a unit value, e.g. "4w" for 4 weeks.
    #   offset - a span string, consisting of a length integer
    #       and a unit value, e.g. "0d" for no offset.

    # Details:
    #   Allowed unit values are "m" for 4 weeks, "w" for weeks,
    #   "d" for days, "H" for hours, "M" for minutes, and "S"
    #   for seconds.

    # Example:
    #   periodicallyRolling((time(as.timeSeries(data(smallcap.ts)))))

    # FUNCTION:

    # Check:
    stopifnot(is(x, "timeDate"))

    # Settings:
    periods = c(4*7*24*3600, 7*24*3600, 24*3600, 3600, 60, 1)
    names(periods) = Names = c("m", "w", "d", "H", "M", "S")
    periodUnit = gsub("[ 0-9]", "", period, perl = TRUE)
    stopifnot(periodUnit %in% Names)
    offsetUnit = gsub("[ 0-9]", "", offset, perl = TRUE)
    stopifnot(offsetUnit %in% Names)
    byUnit = gsub("[ 0-9]", "", by, perl = TRUE)
    stopifnot(byUnit %in% Names)

    # Extract Periods:
    period = as.integer(gsub("[mwdHMS]", "", period, perl = TRUE)) *
        periods[periodUnit]
    offset = as.integer(gsub("[mwdHMS]", "", offset, perl = TRUE)) *
        periods[offsetUnit]
    by = as.integer(gsub("[mwdHMS]", "", by, perl = TRUE)) *
        periods[byUnit]

    # Convert timeDate to GMT-POSIX
    posixGMT = as.POSIXct(
        timeDate(x, zone = x@FinCenter, FinCenter = "GMT"), tz = "GMT")

    # Compute Julian counts (x) and series values (y)
    Origin = as.POSIXct("1970-01-01", tz = "GMT")
    u <- as.integer(difftime(posixGMT, Origin, tz = "GMT", units = "secs"))
    xout = seq(u[1] + offset, u[length(u)], by = by)
    toGMT = Origin + as.integer(xout)
    fromGMT = toGMT - period
    toGMT = toGMT[fromGMT >= posixGMT[1]]
    fromGMT = fromGMT[fromGMT >= posixGMT[1]]
    to = timeDate(toGMT, zone = "GMT", FinCenter = x@FinCenter)
    from = timeDate(fromGMT, zone = "GMT", FinCenter = x@FinCenter)

    # Windows:
    windows = list(from = from, to = to)
    attr(windows, "control") = c(start = start(x), end = end(x))

    # Return Value:
    windows
}


# ---------------------------------------------------------------------------- #
# Roxygen Tags
#' @export
# ---------------------------------------------------------------------------- #
monthlyRolling <-
    function(x, period = "12m", by = "1m")
{
    # A function implemented by Diethelm Wuertz and Yohan Chalabi

    # Description:
    #   Returns start and end dates for monthly periods

    # Arguments:
    #   x - an object of class timeDate
    #   period - a span string, consisting of a length integer
    #       and a unit value, e.g. "12m" for 1 calendar year.
    #   by - a span string, consisting of a length integer
    #       and a unit value, e.g. "1m" for 1 calendar month.

    # Details:
    #   The only allowed allowed unit value is "m" for monthly
    #   periods. Express a quarterly period by "3m", a semester
    #   by "6m", a year by "12m" etc.

    # Example:
    #   monthlyRolling((time(as.timeSeries(data(smallcap.ts)))))

    # FUNCTION:

    # Check:
    stopifnot(is(x, "timeDate"))

    # Get Window Parameter:
    periodLength = as.numeric(substr(period, 1, nchar(period)-1))
    periodUnit = substr(period, nchar(period), nchar(period))
    byLength = as.numeric(substr(by, 1, nchar(by)-1))
    byUnit = substr(by, nchar(by), nchar(by))
    stopifnot(periodUnit == "m")
    stopifnot(byUnit == "m")

    # Make Windows - expand series x to a monthly series
    positions = x
    startPositions = unique(timeFirstDayInMonth(positions))

    # for non monthly data
    # startPositions@Data[1] <- start(x)@Data
    endPositions = unique(timeLastDayInMonth(positions))

    # for non monthly data
    # endPositions@Data[length(endPositions)] <- end(x)@Data
    numberOfPositions = length(startPositions)
    startSeq <- seq(from = 1,
        to = (numberOfPositions-periodLength + 1),
        by = byLength)
    startDates = startPositions[startSeq]
    endSeq <- seq(from = periodLength,
        to = numberOfPositions,
        by = byLength)
    endDates = endPositions[endSeq]

    # Windows:
    windows = list(from = startDates, to = endDates)
    attr(windows, "control") = c(start = start(positions), end = end(positions))

    # Return Value:
    windows
}


################################################################################

