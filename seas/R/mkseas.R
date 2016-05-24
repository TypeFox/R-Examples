"mkseas" <-
function (x, width=11, start.day=1, calendar, year) {
  if (missing(start.day) || is.null(start.day))
    start.day <- 1
  if (missing(calendar) || is.null(calendar)) {
    missing.calendar <- TRUE
    calendar <- NULL
  } else missing.calendar <- FALSE
  if (missing(x)) {
    missing.x <- TRUE
    if (is.null(calendar) || calendar=="julian") {  # have leap years
      if (missing(year))
        stop("a 'year' is needed to determine if this is a leap year")
      year.length <- year.length(year[1], calendar)
    } else year.length <- year.length(1, calendar)
    x <- 1:year.length
    orig <- "x"
  } else {
    missing.x <- FALSE
    orig <- as.character(substitute(x))[[1]]
    if (missing.calendar) {
      if ("date" %in% names(x) && "calendar" %in% names(attributes(x$date)))
        calendar <- attr(x$date, "calendar")
      else if ("calendar" %in% names(attributes(x)))
        calendar <- attr(x, "calendar")
    }
    year.length <- year.length(2000, calendar)
  }
  if (inherits(start.day, c("POSIXct", "Date"))) {
    start.month <- as.integer(format(start.day, "%m"))
    start.std <- as.Date(format(start.day, "2000-%m-%d"))
    start.yday <- as.integer(format(start.std, "%j"))
  } else {  # assumed integer
    start.day <- as.integer(start.day)
    if (start.day < 1 || start.day >= year.length)
      stop(gettextf("%s must >= 1 or < %s",
                    sQuote("start.day"), sQuote("year.length")))
    start.yday <- start.day
    start.std <- as.Date(sprintf("2000-%i", as.integer(start.yday)), "%Y-%j")
    start.month <- as.integer(format(start.std, "%m"))
  }
  if (inherits(x, c("POSIXct", "Date"))) {  # find the date
    date <- x
  } else if (inherits(x, "data.frame")) {
    if (inherits(x$date, c("POSIXct", "Date")))
      date <- x$date
    else {
      cls <- sapply(x, class) %in% c("POSIXct", "Date")
      if (any(cls)) {
        date <- names(x)[cls]
        warning(gettextf("%s column was found in %s",
                         sQuote("date"),
                         sQuote(sprintf("%s$%s", orig, date))))
        date <- x[[date]]
      } else {
        stop(gettextf("could not find a %s colum in %s",
                      sQuote("date"), sQuote(orig)))
      }
    }
  } else if (inherits(x, c("numeric", "integer"))) {  # assumed to be yday
    date <- as.Date(sprintf("%04i-%03i", as.integer(year), as.integer(x)),
                    "%Y-%j")
  } else {
    stop(gettextf("could not find %s in %s", sQuote("date"), sQuote(orig)))
  }
  if (is.character(width)) {
    if (year.length %in% c(365, 366)) {  # 366 days per year
      days <- c(31, 29, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31)
      if (year.length == 365)
        days[2] <- 28  # 365-day calendar
    } else {
      days <- rep(year.length / 12, 12)
    }
    if (width == "mon") {  # short month name from locale
      levels <- months(as.Date(paste(2000, 1:12, 1, sep="-")), TRUE)
      bins <- format(date, "%b")
      start.bin <- start.month
    } else if (width == "month") {  # long month name from locale
      levels <- months(as.Date(paste(2000, 1:12, 1, sep="-")), FALSE)
      bins <- format(date, "%B")
      start.bin <- start.month
    } else if (width == "DJF") {
      days <- c(sum(days[c(1:2, 12)]), sum(days[3:5]),
                sum(days[6:8]), sum(days[9:11]))
      levels <- c("DJF", "MAM", "JJA", "SON")  # these are not locale specific
      bins <- levels[floor((as.integer(format(date, "%m")) %% 12) / 3 ) + 1]
      start.bin <- floor(start.month %% 12 / 3) + 1
    } else if (width == "JF") {
      days <- c(sum(days[1:2]), sum(days[3:4]), sum(days[5:6]),
                sum(days[7:8]), sum(days[9:10]), sum(days[11:12]))
      levels <- c("JF", "MA", "MJ", "JA", "SO", "ND")
      bins <- levels[floor(((as.integer(format(date, "%m")) - 1) %% 12) / 2) + 1]
      start.bin <- floor((start.month - 1) %% 12 / 2) + 1
    } else if (width == "JFM") {
      days <- c(sum(days[1:3]), sum(days[4:6]),
                sum(days[7:9]), sum(days[10:12]))
      levels <- c("JFM", "AMJ", "JAS", "OND")
      bins <- levels[floor(((as.integer(format(date, "%m")) - 1) %% 12) / 3) + 1]
      start.bin <- floor((start.month - 1) %% 12 / 3) + 1
    } else if (width %in% c("zod", "zodiac")) {
      calendar <- NA  # Gregorian
      days <- c(30, 31, 32, 31, 31, 31, 30, 30, 30, 29, 30, 31)
      levels <- c("Aries", "Taurus", "Gemini", "Cancer", "Leo", "Virgo",
                  "Libra", "Scorpio", "Sagittarius", "Capricorn", "Aquarius",
                  "Pisces")
      if (width == "zod")
          levels <- abbreviate(levels, 3)
      # lower dates for each sign in `mmdd' form; starting with Aquarius
      # using http://en.wikipedia.org/wiki/Horoscope#Western_astrological_sign
      zod.b <- c(120, 219, 321, 420, 521, 622, 723, 823, 923, 1023, 1122, 1222)
      md <- as.integer(format(date, "%m%d"))
      # 9 puts Cap in right spot
      bins <- levels[(findInterval(md, zod.b) + 9) %% 12 + 1]
      start.bin <- (findInterval(start.month, zod.b)+9) %% 12 + 1
    } else {  # add other calendar divisions if needed
      stop(gettextf("%s width argument not supported", sQuote(width)))
    }
    n.bins <- length(levels)
    order <- (start.bin:(start.bin - 1 + n.bins) - 1) %% n.bins + 1
    days <- days[order]
    levels <- levels[order]
  } else if (is.numeric(width)) {  # n-day bins
    yday <- as.integer(format(date, "%j"))
    yday[yday > year.length] <- year.length  # trim anything bigger
    ld <- if (is.null(calendar))  # Gregorian
      function(y)(ifelse(y %% 4 == 0 & y %% 100 !=0 | y %% 400 == 0, 0, -1))
    else if (calendar == "julian")
      function(y)(ifelse(y %% 4 == 0, 0, -1))
    else
      function(y)(0)
    if (inherits(x, "data.frame") && "date" %in% names(x)) {
      if (!"year" %in% names(x))
        x$year <- as.integer(format(x$date, "%Y"))
      lds <- ld(x$year)
    } else {
      lds <- integer(year.length)
    }
    if (any(lds!=0) && inherits(start.day, c("POSIXct", "Date"))) {
      if (start.yday > 60) {  # before or after March 1
        yday <- yday-start.yday - lds
      } else {
        if (!is.null(x$date))
          lds <- ld(x$year - 1)
        yday <- yday - start.yday
      }
    } else {
      yday <- yday - start.yday
    }
    bins <- floor(yday %% (year.length + lds) / width + 1)
    nbins <- floor((year.length - 1) / width)  # 0-based
    if ((year.length - 1) / width - nbins < 0.2)
      bins[bins > nbins] <- nbins  # trim last bin
    else
      nbins <- nbins + 1  # now 1-based
    days <- rep(width, nbins)
    days[nbins] <- width + year.length - sum(days)
    levels <- 1:nbins  # 1-based
  }
  bins <- factor(bins, levels)
  attr(bins, "width") <- width
  attr(bins, "start.day") <- start.day
  attr(bins, "calendar") <- calendar
  attr(bins, "bin.lengths") <- days
  bins
}
