### Susumu Tanimura <aruminat@gmail.com>
### Time-stamp: <2015-06-24 17:16:11 umusus>
### This code was inspired by and refferrd to VBA macro at
### http://www.h3.dion.ne.jp/~sakatsu/holiday_logic3.htm

jholiday <- function(year, holiday.names = TRUE){
  if(missing(year)) stop("year was not given")
  stopifnot(is.numeric(year))
  stopifnot(length(year) == 1)
  stopifnot(year > 1948)
  ## We need strings of weekdays for a different locale.
  ## Starts with Sunday
  wdayStrings <- weekdays(seq(as.Date("2013/2/17"), by = "days", length = 7))

  .fixedDate <- function(x, name){
    h <- as.Date(paste(year, x, sep = "-"))
    if(holiday.names) names(h) <- name
    return(c(d, h))
  }
  .unfixedDate <- function(h, name){
    if(holiday.names) names(h) <- name
    return(c(d, h))
  }
  d <- as.Date(vector())                # There are no Date()
  n <- character()
  ## ====== January ==========
  # New Year's Day
  d <- .fixedDate("01-01", "New Year's Day")
  # Coming of Age Day
  if(year >= 2000){
    d <- .unfixedDate(.findDateByWeekday(year, 1, wdayStrings[2], 2),
                      "Coming of Age Day")
  }else{
    d <- .fixedDate("01-15", "Coming of Age Day")
  }

  ## ====== February ==========
  # Foundation Day
  if(year >= 1967){
    d <- .fixedDate("02-11", "Foundation Day")
  }else if(year == 1989){
    d <- .fixedDate("02-24", "State Funeral of the Showa Emperor")
  }

  ## ====== March ==========
  # Vernal Equinox Day
  d <- .unfixedDate(as.Date(paste(year, "03", .Shunbun(year), sep = "-")),
                    "Vernal Equinox Day")

  ## ====== April ==========
  # Showa Day
  if(year >= 2007){
    d <- .fixedDate("04-29", "Showa Day")
  }else if(year >= 1989){
    d <- .fixedDate("04-29", "Greenery Day")
  }else {
    d <- .fixedDate("04-29", "The Emperor's Birthday")
  }
  # Marriage of Crown Prince Akihito
  if(year == 1959){
    d <- .fixedDate("04-10", "Marriage of Crown Prince Akihito")
  }
  ## ====== May ==========
  # Constitution Memorial Day
  d <- .fixedDate("05-03", "Constitution Memorial Day")
  if(year >= 2007){
    d <- .fixedDate("05-04", "Greenery Day")
  }else if(year >= 1988){
    d <- .fixedDate("05-04", "Citizens' Holiday")
  }
  # Children's Day
  d <- .fixedDate("05-05", "Children's Day")
  
  ## ====== June ==========
  if(year == 1993){
    d <- .fixedDate("06-09", "Marriage of Crown Prince Naruhito")
  }

  ## ====== July ==========
  # Marine Day
  if(year >= 2003){
    d <- .unfixedDate(.findDateByWeekday(year, 7, wdayStrings[2], 3),
                      "Marine Day")
  } else if(year >= 1996){
    d <- .fixedDate("07-20", "Marine Day")
  }

  ## ====== August ==========
  # Mountain Day
  if(year >= 2016){
    d <- .fixedDate("08-11", "Mountain Day")
  }
  
  ## ====== September ==========
  # Autumnal Equinox Day
  aed <- as.Date(paste(year, "09", .Shubun(year), sep = "-"))
  d <- .unfixedDate(aed, "Autumnal Equinox Day")
  # Respect-for-the-Aged Day
  if(year >= 2003){
    rad <- .findDateByWeekday(year, 9, wdayStrings[2], 3)
    d <- .unfixedDate(rad, "Respect-for-the-Aged Day")
    if((aed - rad) == 2){
      d <- .unfixedDate(aed - 1, "Citizens' Holiday")
    }
  }else if(year >= 1966){
    d <- .fixedDate("09-15", "Respect-for-the-Aged Day")
  }

  ## ====== October ==========
  # Health and Sports Day
  if(year >= 2000){
    d <- .unfixedDate(.findDateByWeekday(year, 10, wdayStrings[2], 2),
                      "Health and Sports Day")
  }else if(year > 1966){
    d <- .fixedDate("10-10", "Health and Sports Day")
  }
  
  ## ====== November ==========
  # Culture Day
  d <- .fixedDate("11-03", "Culture Day")
  # Labour Thanksgiving Day
  d <- .fixedDate("11-23", "Labour Thanksgiving Day")
  # Official Enthronement Ceremony of Emperor Akihito
  if(year == 1990){
    d <- .fixedDate("11-12", "Official Enthronement Ceremony of Emperor Akihito")
  }

  ## ====== December ==========
  # The Emperor's Birthday
  if(year >= 1989){
    d <- .fixedDate("12-23", "The Emperor's Birthday")
  }

  ## ====== Others ==========
  # Transfer Holiday
  sun <- weekdays(d) == wdayStrings[1]
  if(any(sun)){
    hm <- d[sun] + 1
    if(year >= 2007){
      while(any(i <- hm %in% d)){
        hm[i] <- hm[i] + 1
      }
      if(holiday.names) names(hm) <- rep("Transfer Holiday", length(hm))
      d <- c(d, hm)
    }else if(year >= 1973){
      hm <- hm[! hm %in% d]
      if(holiday.names) names(hm) <- rep("Transfer Holiday", length(hm))
      d <- c(d, hm)
    }
  }
  return(sort(d))
}

is.jholiday <- function(dates){
  if(missing(dates)) stop("dates was not given")
  stopifnot(class(dates) == "Date")
  sapply(dates, function(x){
    y <- as.numeric(format(x, "%Y"))
    h <- jholiday(y, holiday.names = FALSE)
    x %in% h})
}

### ====== Fix me =======
### improve calculation cost
### =====================
.findDateByWeekday <- function(year, month, weekday, ordinal){
  stopifnot(is.numeric(year))
  stopifnot(is.numeric(month))
  stopifnot(is.character(weekday))
  stopifnot(is.numeric(ordinal))
  from <- as.Date(paste(year, month, 1, sep="/"))
  to <- as.Date(paste(year, month + 1, 1, sep="/")) - 1
  d <- seq(from, to, by="days")
  w <- weekdays(d) == weekday
  x <- d[w]
  return(x[ordinal])
}
  
## Calculation formulas in .Shunbun() and .Shubun() were referred to
## Shinkoyomi-benricho, Ephemeris Computation Workshop eds, Koseisha
## Koseikaku: Tokyo, 1991, ISBN: 9784769907008.

.Shunbun <- function(year){
  if(year <= 1947){
    dd = NULL;
  }
  else if(year <= 1979){
    dd = trunc(20.8357 + 0.242194 * (year - 1980) - trunc((year - 1983) / 4))
  }
  else if(year <= 2099){
    dd = trunc(20.8431+ 0.242194 * (year - 1980) - trunc((year - 1980) / 4))
  }
  else if(year <= 2150){
    dd = trunc(21.851 + 0.242194 * (year-1980) - trunc((year - 1980) / 4))
  }
  else {
    dd = NULL
  }
  return(dd)
}

.Shubun <- function(year){
  if(year <= 1947){
    dd = NULL;
  }
  else if(year <= 1979){
    dd = trunc(23.2588 + 0.242194 * (year - 1980) - trunc((year - 1983) / 4))
  }
  else if(year <= 2099){
    dd = trunc(23.2488+ 0.242194 * (year - 1980) - trunc((year - 1980) / 4))
  }
  else if(year <= 2150){
    dd = trunc(24.2488 + 0.242194 * (year-1980) -  trunc((year - 1980) / 4))
  }
  else {
    dd = NULL
  }
  return(dd)
}
