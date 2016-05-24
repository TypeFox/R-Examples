
calendar.effects <- function(x, trading.day = TRUE, easter = 6, 
  leap.year = FALSE, holidays = NULL, easter.date = FALSE)
{
  easterSunday <- function(year)
  #http://code.activestate.com/recipes/576517-calculate-easter-western-given-a-year/
  {
    a <- year %% 19
    b <- floor(year / 100)
    c <- year %% 100
    d <- (19 * a + b - floor(b / 4) - floor((b - floor((b + 8) / 25) + 1) / 3) + 15) %% 30
    e <- (32 + 2 * (b %% 4) + 2 * floor(c / 4) - d - (c %% 4)) %% 7
    f <- d + e - 7 * floor((a + 11 * d + 22 * e) / 451) + 114
    d <- as.Date(paste(year, floor(f / 31), f %% 31 + 1, sep="/"))
    if (weekdays(d) != "Sunday")
      stop("unexpected weekday for Easter day")
    d
  }

  # proportion of days before Easter and Easter Sunday
  # belonging to months March and April

  easter.prop <- function(year, backdays = 6, pmarch = 0.5, papril = 0.5)
  {
    stopifnot(pmarch + papril == 1)
    d1 <- easterSunday(year)
    d2 <- d1 - backdays
    days <- seq(from = d2, to = d1, by = "day")
    #days <- seq(from = d1, to = d2, by = -1)

    list(easter = d1, 
      prop = c(sum(months(days) == "March") / (backdays + 1) - pmarch, 
      sum(months(days) == "April") / (backdays + 1) - papril))
  }

  #stopifnot(is.ts(x))
  s <- frequency(x)
  if (s != s)
    stop("this function is defined only for monthly time series")

  if (!is.null(holidays))
    stopifnot(tsp(x) == tsp(holidays))

  if (!trading.day && !is.null(holidays))
    warning("argument ", sQuote("holidays"), " was ignored since argument ", 
      sQuote("trading.day"), " is FALSE")

  # Easter

  if (easter > 0)
  {
    y0 <- start(x)[1]
    yN <- end(x)[1]
    ereg <- ts(0, start = c(y0, 1), end = c(yN, 12), frequency = 12)
    Edate <- rep(NA, length(yN - y0 + 1))
    id <- matrix(which(cycle(ereg) %in% c(3,4)), ncol = 2, byrow = TRUE)

    i <- 1
    for (y in seq.int(y0, yN))
    {
      tmp <- easter.prop(year = y, backdays = easter)
      Edate[i] <- as.character(tmp$easter)
      ereg[id[1,]] <- tmp$prop
      if (y < yN - 1) {
        id <- id[-1,]
      } else
        id <- rbind(id[-1,])
      i <- i + 1
    }

    ereg <- window(ereg, start = start(x), end = end(x))
  } else
    ereg <- NULL

  # trading day and leap year

  if (trading.day) {
    td <- ts(rep(NA, length(x)))
    tsp(td) <- tsp(x)
  } else 
    td <- NULL

  if (leap.year) {
    ly <- ts(rep(0, length(x)))
    tsp(ly) <- tsp(x)
  } else 
    ly <- NULL

  if (trading.day || leap.year)
  {
    for (i in seq_along(x))
    {
      ti <- time(x)[i]
      fti <- floor(ti)
      si <- round(s * (ti - fti)) + 1
      
      d0 <- as.Date(paste(c(fti, si, 1), collapse = "-"), format = "%Y-%m-%d")

      if (si == s) {
        ysN <- c(fti + 1, 1)
      } else
        ysN <- c(fti, si + 1)

      dN <- as.Date(paste(c(ysN, 1), collapse = "-"), 
        format = "%Y-%m-%d") - 1

      days <- weekdays(seq(from = d0, to = dN, by = "day"), abbreviate = TRUE)

      if (trading.day)
      {
        wd <- sum(days %in% c("Mon", "Tue", "Wed", "Thu", "Fri"))
        nwd <- sum(days %in% c("Sat", "Sun"))

        if (!is.null(holidays))
        {
          wd <- wd - holidays
          nwd <- nwd + holidays
        }

        # number of working days minus number of non-working days x 5/2

        td[i] <- wd - (5/2) * nwd
      }

      if (leap.year)
      {
        if (si == 2)
        {
          ndays <- length(days)
          if (ndays == 28) {
            ly[i] <- -0.25
          } else 
          if (ndays == 29) {
            ly[i] <- 0.75
          } else
            stop("unexpected number of days in February")
        } #else
          #ly[i] <- 0
      }
    }
  }

  m <- cbind(trading.day = td, leap.year = ly, easter = ereg)
  if (is.null(dim(m))) {
    m <- ts(matrix(m))
    tsp(m) <- tsp(x)
  }

  colnames(m) <- c("trading-day", "leap-year", "Easter")[#
    c(trading.day, leap.year, easter > 0)]

  if (easter.date) {
    return(list(effects = m, easter = Edate))
  } else
    return(m)
}
