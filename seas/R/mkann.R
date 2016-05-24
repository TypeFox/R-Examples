"mkann" <-
function(x, start.day, calendar) {
  orig <- as.character(substitute(x))[[1]]
  if (inherits(x, c("POSIXct", "Date"))) {  # find the date
    x <- data.frame(date=x)
  } else if (inherits(x, "data.frame")) {
    if (!inherits(x$date, c("POSIXct", "Date"))) {
      cls <- sapply(x, class) %in% c("POSIXct", "Date")
      if (any(cls)) {
        date <- names(x)[cls]
        warning(gettextf("%s column was found in %s",
                         sQuote("date"),
                         sQuote(sprintf("%s$%s", orig, date))))
        x$date <- x[[date]]
      } else {
        stop(gettextf("could not find a %s colum in %s",
                      sQuote("date"), sQuote(orig)))
      }
    }
  } else {
    stop(gettextf("could not find %s in %s", sQuote("date"), sQuote(orig)))
  }
  if (missing(start.day) || is.null(start.day))
    start.day <- c(attr(x, "start.day"), 1)[[1]]
  if (missing(calendar) || is.null(calendar)) {
    calendar <- if ("date" %in% names(x) && "calendar" %in% names(attributes(x$date)))
      attr(x$date, "calendar")
    else if ("calendar" %in% names(attributes(x)))
      attr(x, "calendar")
    else
      NULL
  }
  year.length <- year.length(2000, calendar)
  year <- as.integer(format(x$date, "%Y"))
  year.range <- range(year)
  days <- if (year.length %in% c(365, 366))
    c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
  else
    rep(year.length / 12, 12)
  if (year.length == 365)
    days[2] <- 28  # non-Gregorian calender
  if (inherits(start.day, c("POSIXct", "Date"))) {
    start.std <- as.Date(format(start.day, "2000-%m-%d"))
    start.yday <- as.integer(format(start.std, "%j"))
  } else {
    start.yday <- start.day
  }
  ld <- if (is.null(calendar))  # Gregorian)
    function(y)(ifelse(y %% 4 == 0 & y %% 100 != 0 | y %% 400 == 0, 0, -1))
  else if (calendar == "julian")
    function(y)(ifelse(y %% 4==0, 0, -1))
  else
    function(y)(integer(length(y)))
  if (start.yday != 1) {  # not algned to years
    yday <- as.integer(format(x$date, "%j"))
    yday[yday > year.length] <- year.length
    if (inherits(start.day, c("POSIXct", "Date"))) {
      lds <- ld(year)
      if (any(lds!=0)) {
        if (start.yday > 60) {  # adjust after March 1st
          yday <- yday - start.yday - lds
        } else {
          yday <- yday - start.yday
        }
      }
    } else {
      if (start.day < 1 || start.day >= year.length)
        stop(gettextf("%s must >= 1 or < %s",
                      sQuote("start.day"), sQuote("year.length")))
      yday <- yday - start.day
    }

    flr <- floor(yday / (year.length + lds))
    year <- year + flr
    if (year.length == 366 && start.yday > 60)
      year.lengths <- year.length + ld(unique(year + 1))
    else
      year.lengths <- year.length + ld(unique(year))
    year <- paste(year, year + 1, sep="_")
  } else {
    year.lengths <- year.length + ld(unique(year))
  }
  x <- factor(year)
  attr(x, "start.day") <- start.day
  attr(x, "year.range") <- year.range
  attr(x, "calendar") <- calendar
  attr(x, "year.lengths") <- year.lengths
  x
}
