.SetTimeZero <- function(time0, y) {
  if (is.null(time0)) {
    if (is.null(y)) {
      stop("You must supply time0 if y is missing.")
    }
    if (!inherits(y, "zoo")) {
      ## Note:  an xts object inherits from zoo.
      stop("You must supply 'time0' if y is not a zoo or xts object.")
    }
    times <- index(as.xts(y))
    tryCatch(time0 <- as.POSIXct(times)[1],
             error = simpleError(
               "The index of y could not be converted to POSIXt."))
  }
  stopifnot(inherits(time0, "POSIXt"))
  return(as.POSIXlt(time0))
}

.ValidateHolidaySigmaPrior <- function(sigma.prior, sdy) {
  ## The prior distribution says that sigma is small, and can be no
  ## larger than the sample standard deviation of the time series
  ## being modeled.
  if (is.null(sigma.prior)) {
    return(SdPrior(.01 * sdy, upper.limit = sdy))
  }
  stopifnot(inherits(sigma.prior, "SdPrior"))
  return(sigma.prior)
}

.ValidateHolidayInitialStatePrior <- function(initial.state.prior, sdy) {
  if (is.null(initial.state.prior)) {
    return(NormalPrior(0, sdy))
  }
  stopifnot(inherits(initial.state.prior, "NormalPrior"))
  return(initial.state.prior)
}

.BaseHolidaySpecification <- function(holiday.name,
                                      time0,
                                      y,
                                      sdy,
                                      sigma.prior,
                                      initial.state.prior,
                                      days.before,
                                      days.after) {
  time0 <- .SetTimeZero(time0, y)
  initial.state.prior <-
    .ValidateHolidayInitialStatePrior(initial.state.prior, sdy)
  sigma.prior <- .ValidateHolidaySigmaPrior(sigma.prior, sdy)
  stopifnot(is.character(holiday.name) && length(holiday.name) == 1)
  stopifnot(is.numeric(days.before) &&
            length(days.before) == 1 &&
            days.before >= 0)
  stopifnot(is.numeric(days.after) &&
            length(days.after) == 1 &&
            days.after >= 0)

  holiday.spec <- list(name = holiday.name,
                       sigma.prior = sigma.prior,
                       sdy = sdy,
                       sigma.prior = sigma.prior,
                       initial.state.prior = initial.state.prior,
                       days.before = days.before,
                       days.after = days.after,
                       time0.month = as.integer(1 + time0$mon),
                       time0.day = as.integer(time0$mday),
                       time0.year = as.integer(1900 + time0$year),
                       size = 1 + days.before + days.after)
  return(holiday.spec)
}

AddFixedDateHoliday <- function(state.specification = NULL,
                                holiday.name,
                                month,
                                day,
                                y,
                                sigma.prior = NULL,
                                initial.state.prior = NULL,
                                sdy = sd(as.numeric(y), na.rm = TRUE),
                                time0 = NULL,
                                days.before = 1,
                                days.after = 1) {
  ## Adds a fixed date holiday to the state specification.  A fixed
  ## date holiday is a holiday that occurs on the same date each year
  ## (e.g. Christmas is December 25).
  ##
  ## Args:
  ##   state.specification: A list of state components.
  ##   holiday.name:  The name of the holiday (a string).
  ##   month: The name of the month in which the holiday occurs.  The
  ##     first letter should be capitalized.  Unambiguous partial
  ##     matches are allowed.
  ##   day:  The day of the month on which the holiday occurs.
  ##   y: The time series to be modeled, preferably as a zoo or xts
  ##     object.  Must be daily data.
  ##   sigma.prior: Prior distribution on the standard deviation of
  ##     class SdPrior.
  ##   initial.state.prior: An object created by NormalPrior.
  ##     The prior distribution on the values of the initial state
  ##     (i.e. the state of the first observation).
  ##   sdy: The standard deviation of y.  This will be ignored if y is
  ##     provided, or if both sigma.prior and initial.state.prior are
  ##     supplied directly.
  ##   time0: An object convertible to class POSIXt (such as a Date
  ##     or character string) giving the date of the first day in the
  ##     training data y.
  ##   days.before: The number of days the influence of the named
  ##     holidays extends prior to the actual holiday.
  ##   days.after: The number of days the influence of the named
  ##     holidays extends after the actual holiday.
  ##
  ## Returns:
  ##   The requested holidays are appended to state.specification, and
  ##   the modified state.specification is returned.
  if (is.null(state.specification)) {
    state.specification <- list()
  }
  stopifnot(is.list(state.specification))
  if (missing(y)) y <- NULL
  holiday.spec <- .BaseHolidaySpecification(holiday.name,
                                            time0,
                                            y,
                                            sdy,
                                            sigma.prior,
                                            initial.state.prior,
                                            days.before,
                                            days.after)
  day <- as.integer(day)
  stopifnot(day >= 1 && day <= 31)
  holiday.spec$day <- day
  month <- match.arg(month, month.name)
  holiday.spec$month <- month
  class(holiday.spec) <- c("FixedDateHoliday", "Holiday", "StateModel")
  state.specification[[length(state.specification) + 1]] <- holiday.spec
  return(state.specification)
}

AddNthWeekdayInMonthHoliday <- function(state.specification = NULL,
                                        holiday.name,
                                        month,
                                        day.of.week,
                                        which.week,
                                        y,
                                        sigma.prior = NULL,
                                        initial.state.prior = NULL,
                                        sdy = sd(as.numeric(y), na.rm = TRUE),
                                        time0 = NULL,
                                        days.before = 1,
                                        days.after = 1) {
  ## A holiday that occurs on a fixed occurrence of a particular week
  ## day in the same month every year.  E.g. the third Monday in
  ## January.
  ## Args:
  ##   holiday.name:  The name of the holiday.
  ##   month: The month in which the holiday takes place, as a string.
  ##     Abbreviations are okay.  Capitalize the first letter.
  ##   day.of.week: The name of the day of the week when the holiday
  ##     takes place.  Abbreviations are okay.  Capitalize the first
  ##     letter.
  ##   which.week: An integer specifying the week in which the holiday
  ##     takes place.  If which.week <= 0 then holiday occurs on the
  ##     last specified day.of.week in the month.
  ##   sigma.prior: Prior distribution on the standard deviation of
  ##     class SdPrior.
  ##   initial.state.prior: An object created by NormalPrior.
  ##     The prior distribution on the values of the initial state
  ##     (i.e. the state of the first observation).
  ##   sdy: The standard deviation of y.  This will be ignored if y is
  ##     provided, or if both sigma.prior and initial.state.prior are
  ##     supplied directly.
  ##   time0: An object convertible to class POSIXt (such as a Date
  ##     or character string) giving the date of the first day in the
  ##     training data y.
  ##   days.before: The number of days the influence of the named
  ##     holidays extends prior to the actual holiday.
  ##   days.after: The number of days the influence of the named
  ##     holidays extends after the actual holiday.
  ##
  ## Example:
  ##   To specify the third Monday in January,
  ##   AddNthWeekdayInMonthHoliday(list(), "NationalSteveDay", "January",
  ##       "Monday", 3)
  ##
  ## Returns:
  ##   An object of class NthWeekdayInMonthHoliday.
  if (is.null(state.specification)) {
    state.specification <- list()
  }
  stopifnot(is.list(state.specification))
  if (missing(y)) y <- NULL
  holiday.spec <- .BaseHolidaySpecification(holiday.name,
                                            time0,
                                            y,
                                            sdy,
                                            sigma.prior,
                                            initial.state.prior,
                                            days.before,
                                            days.after)
  month <- match.arg(month, month.name)
  holiday.spec$month <- month
  day.of.week.names <- c("Sunday", "Monday", "Tuesday", "Wednesday",
                         "Thursday", "Friday", "Saturday")
  day.of.week <- match.arg(day.of.week, day.of.week.names)
  holiday.spec$day.of.week <- day.of.week
  which.week <- as.integer(which.week)
  stopifnot(which.week <= 5)
  holiday.spec$which.week <- which.week
  class(holiday.spec) <-
    c("NthWeekdayInMonthHoliday", "Holiday", "StateModel")
  if (which.week <= 0) {
    class(holiday.spec) <- c("LastWeekdayInMonthHoliday", class(holiday.spec))
  }
  state.specification[[length(state.specification) + 1]] <- holiday.spec
  return(state.specification)
}

AddLastWeekdayInMonthHoliday <- function(state.specification = NULL,
                                         holiday.name,
                                         month,
                                         day.of.week,
                                         y,
                                         sigma.prior = NULL,
                                         initial.state.prior = NULL,
                                         sdy = sd(as.numeric(y), na.rm = TRUE),
                                         time0 = NULL,
                                         days.before = 1,
                                         days.after = 1) {
  ## Syntactic sugar for an NthWeekdayInMonthHoliday with N set to -1.
  return(AddNthWeekdayInMonthHoliday(state.specification,
                                     holiday.name,
                                     month,
                                     day.of.week,
                                     which.week = -1,
                                     y,
                                     sigma.prior,
                                     initial.state.prior,
                                     sdy,
                                     time0,
                                     days.before,
                                     days.after))
}


NamedHolidays <- function(except = NULL) {
  ## Returns a list of the named holidays that can be used with
  ## AddNamedHolidays.
  ## Args:
  ##   except: A character vector giving the names of holidays to
  ##     omit.  Matching is done using grep, so any partial match to
  ##     any holiday name will be excluded.
  ##
  ## Returns:
  ##   A character vector giving the full names of the holidays.
  holidays <- c("NewYearsDay",
                "MartinLutherKingDay",
                "SuperBowlSunday",
                "PresidentsDay",
                "ValentinesDay",
                "SaintPatricksDay",
                "USDaylightSavingsTimeBegins",
                "USDaylightSavingsTimeEnds",
                "EasterSunday",
                "USMothersDay",
                "IndependenceDay",
                "LaborDay",
                "ColumbusDay",
                "Halloween",
                "Thanksgiving",
                "MemorialDay",
                "VeteransDay",
                "Christmas")
  if (!is.null(except)) {
    ## Note I tried pmatch, match, and grep here.  None gave the
    ## semantics I'm looking for.
    positions <- sapply(except, function(x) grep(x, holidays))
    if (is.list(positions)) {
      positions <- unlist(positions)
    }
    positions <- unique(positions)
    if (any(is.na(positions))) {
      warning("Argument '",
              except[is.na(positions)],
              "' was not matched in NamedHolidays")
      positions <- positions[!is.na(positions)]
    }
    if (length(positions) > 0) {
      holidays <- holidays[-positions]
    }
  }
  return(holidays)
}

AddNamedHolidays <- function(state.specification = NULL,
                             named.holidays = NamedHolidays(),
                             y,
                             sigma.prior = NULL,
                             initial.state.prior = NULL,
                             sdy = sd(as.numeric(y), na.rm = TRUE),
                             time0 = NULL,
                             days.before = 1,
                             days.after = 1) {
  ## Args:
  ##   state.specification: A list of state components.  If omitted,
  ##     an empty list is assumed.
  ##   named.holiday: A character vector giving the names of one or
  ##     more holidays to include.
  ##   y: A numeric vector containing the time series to be modeled.
  ##      The time series must contain daily data.
  ##   sigma.prior: An object created by SdPrior.  This is the
  ##     prior distribution on the standard deviation of the random walk
  ##     increments for each day in a holiday's window of influence.
  ##   initial.state.prior: An object created by NormalPrior.
  ##     The prior distribution on the values of the initial state
  ##     (i.e. the state of the first observation).
  ##   sdy: The standard deviation of y.  This will be ignored if y is
  ##     provided, or if both sigma.prior and initial.state.prior are
  ##     supplied directly.
  ##   time0: An object convertible to class POSIXt (such as a Date
  ##     or character string) giving the date of the first day in the
  ##     training data y.
  ##   days.before: The number of days the influence of the named
  ##     holidays extends prior to the actual holiday.
  ##   days.after: The number of days the influence of the named
  ##     holidays extends after the actual holiday.
  ##
  ## Returns:
  ##   The requested holidays are appended to state.specification, and
  ##   the modified state.specification is returned.
  if (is.null(state.specification)) {
    state.specification <- list()
  }
  stopifnot(is.list(state.specification))
  if (missing(y)) y <- NULL

  named.holidays <- match.arg(named.holidays, several.ok = TRUE)
  if (length(named.holidays) > 1) {
    for (holiday in named.holidays) {
      state.specification <-
        AddNamedHolidays(state.specification,
                         holiday,
                         y = y,
                         sigma.prior = sigma.prior,
                         initial.state.prior = initial.state.prior,
                         sdy = sdy,
                         time0 = time0,
                         days.before = days.before,
                         days.after = days.after)
    }
  } else if (length(named.holidays) == 1) {
    holiday.spec <- .BaseHolidaySpecification(named.holidays,
                                              time0,
                                              y,
                                              sdy,
                                              sigma.prior,
                                              initial.state.prior,
                                              days.before,
                                              days.after)
    class(holiday.spec) <- c("NamedHoliday", "Holiday", "StateModel")
    state.specification[[length(state.specification) + 1]] <- holiday.spec
  } else {
    stop("Zero-length argument 'named.holidays'")
  }
  return(state.specification)
}
