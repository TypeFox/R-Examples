# Copyright 2011 Google Inc. All Rights Reserved.
# Author: stevescott@google.com (Steve Scott)

## This file contains utility functions for handling date related
## tasks that come up when modeling mixed frequency data.

## TODO(stevescott): Move the date handling functions to a googledates
## package.

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetFractionOfDaysInInitialMonth <- function(week.ending) {
  ## Compute the fraction of days in a week that occur in the month
  ## containing the first day of the week.
  ## Args:
  ##   week.ending: a Date object giving the last date in a week.
  ## Returns:
  ##   A numeric vector of the same length as 'week.ending'.  Each
  ##   entry gives the fraction of days in the week that occur in the
  ##   month containing the start of the week (i.e the date 6 days
  ##   before).
  stopifnot(inherits(week.ending, "Date"))
  begin <- as.POSIXlt(week.ending - 6);
  end <- as.POSIXlt(week.ending)

  fraction.in.initial.month <- rep(1, length(week.ending))
  same.month <- begin$mon == end$mon
  if (any(!same.month)) {
    days <- end$mday[!same.month]
    fraction.in.initial.month[!same.month] <- 1 - days / 7
  }
  return(fraction.in.initial.month)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
GetFractionOfDaysInInitialQuarter <- function(week.ending) {
  ## Compute the fraction of the number of days in a week that occur
  ## in the quarter containing the first day of the week.
  ## Args:
  ##   week.ending:  a Date object giving the last date in a week.
  ## Returns:
  ##   A numeric vector of the same length as 'week.ending'.  Each
  ##   entry gives the fraction of days in the week that occur in the
  ##   quarter containing the start of the week (i.e the date 6 days
  ##   before).
  stopifnot(inherits(week.ending, "Date"))
  begin <- as.POSIXlt(week.ending - 6)
  end <- as.POSIXlt(week.ending)

  same.quarter <- Quarter(begin) == Quarter(end);
  fraction.in.initial.quarter <- rep(1, length(week.ending))
  if (any(!same.quarter)) {
    days <- end$mday[!same.quarter]
    fraction.in.initial.quarter[!same.quarter] <- 1 - days / 7
  }
  return(fraction.in.initial.quarter)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
MatchWeekToMonth <- function(week.ending, origin.month) {
  ## Args:
  ##   week.ending: A vector of class 'Date'.  Each entry is the
  ##     last day in a week.
  ##   origin.month: A scalar vector of class 'Date' that occurs in "month 1".
  ## Returns:
  ##   The index of the month matching the month containing the first
  ##   day in week.ending.  The origin is month 1.  It is the
  ##   caller's responsibility to ensure that these indices correspond
  ##   to legal values in a particular vector of months.
  first.day <- week.ending - 6

  stopifnot(inherits(origin.month, "Date"))
  stopifnot(inherits(week.ending, "Date"))

  if (length(origin.month) > 1) {
    warning("Multiple values supplied for 'origin.month' in ",
            "'MatchWeekToMonth'.  Taking the fist.")
    origin.month <- origin.month[1]
  }

  ## MonthDistance will return the zero-offset distance from the
  ## origin.  Add +1 to get R's unit offset distance.
  return(1 + MonthDistance(first.day, origin.month))
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
WeekEndsMonth <- function(week.ending) {
  ## Indicates which weeks contain the end of a month.
  ## Args:
  ##   week.ending:  A vector of class Date.
  ## Returns:
  ##   A logical vector that is TRUE when the week ending on
  ##   'week.ending' contains the last day in a month, and FALSE
  ##   otherwise.
  first.day <- week.ending - 6
  start.month <- months(first.day)
  end.month <- months(week.ending+1)
  return(end.month != start.month)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
Quarter <- function(date) {
  ## Returns the year and quarter of the year that 'date' belongs to.
  ## The answer is a number of years since 1900, with decimal quarters
  ## in c(.00, .25, .50, .75).
  date <- as.POSIXlt(date)
  quarter <- floor(date$mon / 3)
  year <- date$year
  return(year + quarter / 4)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
WeekEndsQuarter <- function(week.ending) {
  ## Indicates which weeks contain the end of a quarter.
  ## Args:
  ##   week.ending:  A vector of class Date.
  ## Returns:
  ##   A logical vector that is TRUE when the week ending on
  ##   'week.ending' contains the last day in a quarter, and FALSE
  ##   otherwise.  A quarter is defined as ending on the last day of
  ##   March, June, September, or December.
  first.day <- week.ending - 6
  start.quarter <- Quarter(first.day)
  end.quarter <- Quarter(week.ending + 1)
  return(end.quarter != start.quarter)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
MonthDistance <- function(dates, origin) {
  ## Args:
  ##   dates:  A vector of class Date.
  ##   origin:  A singleton vector of class Date giving the origin.
  ## Returns:
  ##   A numeric vector giving the integer number of months between
  ##   'dates' and 'origin'.  Days are ignored, so that if 'dates' is
  ##   in the same month as 'origin' the distance is zero.  The
  ##   distance is signed, so the distance to the preceding month is
  ##   -1, and the same month the preceding year is -12.
  stopifnot(length(origin) == 1)
  origin <- as.POSIXlt(origin)
  dates <- as.POSIXlt(dates)
  ans <- 12 * (dates$year - origin$year) + (dates$mon - origin$mon)
  return(ans)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
ExtendTime <- function(dates, number.of.periods, dt = NULL) {
  ## Extends a sequence of dates to a given length.
  ## Args:
  ##   dates:  A vector of class 'Date'.
  ##   number.of.periods: An integer giving the desired length of the
  ##     output.
  ##   dt: A character string expressing the time lag between elements
  ##     of 'dates'.  Can be "daily", "weekly", "monthly",
  ##     "quarterly", "yearly", or "other".
  ## Returns:
  ##   If number.of.periods is longer than length(dates) then dates is
  ##   padded with extra dates at the end.  Otherwise, 'dates' is
  ##   returned unmodified.
  extra <- number.of.periods - length(dates)
  if (extra <= 0) return(dates)
  if (is.null(dt)) {
    dt <- EstimateTimeScale(dates)
  }
  last.date <- as.Date(tail(dates, 1))

  ## Note that seq.Date behaves strangely when dates are specified by
  ## the last date in a month:
  ##     This works fine:
  ## > seq.Date(from = as.Date("2008-01-20") , by = "month", length.out = 5)
  ## [1] "2008-01-20" "2008-02-20" "2008-03-20" "2008-04-20" "2008-05-20"
  ##     This is wrong:
  ##  > seq.Date(from = as.Date("2008-01-31") , by = "month", length.out = 5)
  ## [1] "2008-01-31" "2008-03-02" "2008-03-31" "2008-05-01" "2008-05-31"
  ##
  ## Thus, when working with monthly dates, we add 1 to the 'from'
  ## argument, and subtract 1 from the resulting sequence:
  ##    This even gets leap year right:
  ## > seq.Date(from = as.Date("2008-01-31") + 1, by = "month", length = 5) - 1
  ## [1] "2008-01-31" "2008-02-29" "2008-03-31" "2008-04-30" "2008-05-31"
  if (dt == "daily") {
    pad <- seq.Date(last.date, by = "day", length.out = extra + 1)
  } else if (dt == "weekly") {
    pad <- seq.Date(last.date, by = "week", length.out = extra + 1)
  } else if (dt == "monthly") {
    pad <- seq.Date(as.Date(last.date) + 1, by = "month",
                    length.out = extra + 1) - 1
  } else if (dt == "quarterly") {
    pad <- seq.Date(last.date + 1, by = "3 months", length.out = extra + 1) - 1
  } else if (dt == "yearly") {
    pad <- seq.Date(last.date, by = "year", length.out = extra + 1)
  } else {
    dt <- mean(diff(dates))
    pad <- last.date + (0:extra) * dt
  }
  dates <- c(head(dates, -1), pad)
  return(dates)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
LastDayInMonth <- function(dates) {
  ## Find the last day in the month containing each date.
  ## Args:
  ##   dates:  A vector object convertible to POSIXlt.  Dates and text
  ## Returns:
  ##   A vector of class 'Date' where each entry is the last day in
  ##   the month containing each entry in 'dates'.
  dates <- as.POSIXlt(dates)
  dates$mday <- 1
  dates$mon <- dates$mon + 1
  next.year <- dates$mon > 11
  if (any(next.year)) {
    dates$mon[next.year] <- 0
    dates$year[next.year] <- dates$year[next.year] + 1
  }
  dates <- as.Date(dates) - 1
  return(dates)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
.Near <- function(x, y, tolerance = 1e-7) {
  ## Rethers whether x is within 'tolerance' of y.  R probably has
  ## something that does this already, but all.equal tried to be too
  ## cute.
  return(abs(x - y) <= tolerance)
}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
EstimateTimeScale <- function(dates) {
  ## Determine the likely time interval between successive dates.
  ## Args:
  ##   dates:  A vector of class 'Date'.
  ## Returns:
  ##   A character string indicating whether the data are "daily",
  ##   "weekly", "monthly", "quarterly", "yearly", or "other".
  dates <- as.Date(dates)
  dt <- as.numeric(diff(dates))
  ndt <- length(dt)
  if (all(.Near(dt, 1.0))) {
    return("daily")
  } else if (all(.Near(dt, 7))) {
    return("weekly")
  } else if (all(.Near(dt, 365, tolerance = 1))) {
    return("yearly")
  } else if (all(.Near(dt, 30, tolerance = 2))) {
    return("monthly")
  } else if (all(.Near(dt, 90, tolerance = 3))) {
    return("quarterly")
  }
  return("other")
}

## MatchDates <- function(fine.dates, coarse.dates) {
##   if (max(fine.dates) > max(coarse.dates)) {
##     extra.dates <- fine.dates[fine.dates > max(coarse.dates)]
##     coarse.dates <- c(coarse.dates, extra.dates)
##   }
##   ans <- match(fine.dates, coarse.dates)
##   if (any(is.na(ans))) {
##     stop("Some dates not matched in MatchDates")
##   }
##   return(ans)
## }

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
AggregateWeeksToMonths <- function(weekly.series,
                                   membership.fraction = NULL,
                                   trim.left = TRUE,
                                   trim.right = NULL) {
  ## A convenience function for aggregating weekly observations to
  ## monthly observations.
  ## Args:
  ##   weekly.series: A zoo time series indexed by the last date in
  ##     each week.  The index must be convertible to 'Date'.
  ##     Multiple time series can be aggregated simultaneously
  ##   membership.fraction: An optional numeric vector giving the
  ##     proportion of a week's measurement attributable to the month
  ##     containing the week's first day.  If NULL then weeks will be
  ##     apportioned to months based on the fraction of the week's
  ##     days that occurs in each month.
  ##   trim.left: Logical indicating whether the first observation in
  ##     the coarse aggregation should be removed.
  ##   trim.right: Logical indicating whether the final observation in
  ##     the coarse aggregate should be removed.  The default behavior
  ##     is to trim unless right endpoint corresponds exactly to the
  ##     end of a coarse interval.
  ## Returns:
  ##   A zoo matrix (if weekly.series was a matrix) or vector
  ##   (otherwise) containing the aggregated time series.  The indices
  ##   of the return value are dates (class Date) set at the last day
  ##   of the month the observation is measuring.
  stopifnot(is.zoo(weekly.series))
  dates <- as.Date(index(weekly.series))

  if (is.null(membership.fraction)) {
    membership.fraction <- GetFractionOfDaysInInitialMonth(index(weekly.series))
  }
  contains.end <- WeekEndsMonth(dates)

  ans <- AggregateTimeSeries(weekly.series,
                             contains.end,
                             membership.fraction,
                             trim.left,
                             trim.right,
                             byrow = TRUE)
  if (is.matrix(ans)) {
    number.of.months <- nrow(ans)
    colnames(ans) <- colnames(weekly.series)
  } else {
    number.of.months <- length(ans)
  }

  # Add date labels to the matrix or vector represented by ans.  The
  # right place to start depends on whether or not the left end point
  # was trimmed, and whether the first week in the series overlapped
  # with a preceding month.  Remember that weeks and months are
  # labelled by their last day.
  initial.date <- as.POSIXlt(dates[1])
  if (trim.left) {
    if (initial.date$mday < 7) {
      ## First week overlaps an earlier month, so the first month is
      ## the month preceding dates[1].  This gets trimmed, which
      ## leaves the month in dates[1] as first.
      initial.date <- LastDayInMonth(dates[1])
    } else {
      ## First week is in the interior of a month that is going to get
      ## trimmed, so initial.date is set to the following month.
      initial.date <- LastDayInMonth(LastDayInMonth(dates[1]) + 1)
    }
  } else {
    if (initial.date$mday < 7) {
      ## The first week overlaps an earlier month, so the first month
      ## label is the month preceding the month containing dates[1].
      initial.date$mday <- 1
      initial.date <- as.Date(initial.date) - 1
    } else {
      ## The first week is entirely in the interior of the month, so
      ## the first month is the same as the month containing dates[1].
      initial.date <- LastDayInMonth(initial.date)
    }
  }

  monthly.dates <- seq.Date(initial.date + 1, by = "month",
                           length.out = number.of.months) - 1
  ans <- zoo(ans, monthly.dates)

}

##=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
AggregateTimeSeries <- function(fine.series,
                                contains.end,
                                membership.fraction,
                                trim.left = any(membership.fraction < 1),
                                trim.right = NULL,
                                byrow = TRUE) {
  ## Aggregate measurements from a fine scaled time series into a
  ## coarse time series.  This is similar to functions from the xts
  ## package, but it can handle weeks -> months aggregation.
  ## Args:
  ##   fine.series: A numeric vector or matrix giving the fine scale
  ##     time series to be aggregated.
  ##   contains.end: A logical vector of the same length as
  ##     fine.series, indicating whether each fine time interval
  ##     contains the end of a coarse time interval.
  ##   membership.fraction: The fraction of each time interval's
  ##     observation attributable to the coarse interval containing
  ##     the fine interval's first day.
  ##   trim.left: Logical indicating whether the first observation in
  ##     the coarse aggregation should be removed.
  ##   trim.right: Logical indicating whether the final observation in
  ##     the coarse aggregate should be removed.  The default behavior
  ##     is to trim unless right endpoint corresponds exactly to the
  ##     end of a coarse interval.
  ##   byrow: Logical.  If fine.series is a matrix, this argument
  ##     indicates whether rows (TRUE) or columns (FALSE) correspond
  ##     to time points.
  ## Returns:
  ##   A matrix (if fine.series is a matrix) or vector (otherwise)
  ##   containing the aggregated time series.
  ##
  ##   Note that unless fine.series happens to coincide with the
  ##   exact beginning or end of a coarse time interval, the left and
  ##   right end points of the resulting aggregation may not contain
  ##   full aggregates.  Use the arguments trim.left and trim.right
  ##   to remove undesired partial aggregates.
  stopifnot(is.numeric(fine.series) || is.matrix(fine.series))
  stopifnot(is.logical(contains.end) && is.numeric(membership.fraction))
  stopifnot(max(membership.fraction) <= 1.0 && min(membership.fraction) > 0)

  transposed <- FALSE
  if (is.matrix(fine.series) && byrow) {
    transposed <- TRUE
    fine.series <- t(fine.series)
  }

  if (is.matrix(fine.series)) {
    time.dimension <- ncol(fine.series)
  } else {
    time.dimension <- length(fine.series)
  }

  stopifnot(time.dimension == length(contains.end) &&
            time.dimension == length(membership.fraction))

  if (is.null(trim.right)) {
    no.remainder <- (tail(contains.end, 1) &&
                     tail(membership.fraction, 1) >= .9999)
    trim.right <- !no.remainder
  }

  aggregate <- .Call("bsts_aggregate_time_series_",
                     fine.series,
                     contains.end,
                     membership.fraction,
                     PACKAGE = "bsts")
  if (is.matrix(aggregate)) {
    if (trim.left) {
      aggregate <- aggregate[, -1, drop = FALSE]
    }
    if (trim.right) {
      aggregate <- aggregate[, -ncol(aggregate), drop = FALSE]
    }
  } else {
    if (trim.left) {
      aggregate <- aggregate[-1]
    }
    if (trim.right) {
      aggregate <- head(aggregate, -1)
    }
  }

  if (transposed) {
    aggregate <- t(aggregate)
  }

  return(aggregate)
}
