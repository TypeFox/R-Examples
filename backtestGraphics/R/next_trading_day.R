#' Next trading day 
#'
#' This function gets the next trading day as defined by is.trading.day()
#'
#' @param d Date relative to which to get the next trading day. 
#' 
#' @return Date The next trading day of date d.

next_trading_day <- function(d) {
  stopifnot(inherits(d, "Date"))
  week <- seq(d + 1, d + 7, by = 1)
  as.Date(week[is_trading_day(week)][1])
}

## Check if the date d is a trading day 

is_trading_day <- function(d) {
  stopifnot(inherits(d, "Date"))
  
  ! weekdays(d) %in% c("Saturday", "Sunday")
  
}