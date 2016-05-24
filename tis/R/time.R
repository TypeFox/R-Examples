is.time <- function(x) all(between(x, 1599, 2500))

time.jul <- function(x, ...){
  y <- jul2ymd(x) %/% 10000
  day <- unclass(x - ymd2jul(10000*y + 101))
  daysInYear <- 365 + isLeapYear(y)
  return(y + day/daysInYear)
}

time.ssDate <- function(x, offset = 1, ...) time(jul(x, offset))

time.ti <- function(x, offset = 1, ...) time(jul(x, offset = offset))

time.POSIXlt <- function(x, ...){
  year <- x$year + 1900
  daysInYear <- 365 + isLeapYear(year)
  year + (x$yday + x$hour/24 + x$min/1440 + x$sec/86400)/daysInYear
}

time.POSIXct <- function(x, ...) time(as.POSIXlt(x))

time2jul <- function(time){
  rawTime <- unclass(time)
  zeros <- numeric(length(rawTime))
  y <- floor(rawTime + 1e-8)
  secondsInYear <- 86400*(zeros + 365 + isLeapYear(y))
  d <- floor(secondsInYear*(rawTime - y) + 0.5)/86400
  return(ymd2jul(10000*y + 0101) + d)
}
