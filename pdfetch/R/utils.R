# Return end of annual date
year_end <- function(date) {
  as.Date(ISOdate(year(date),12,31))
}

# Return end of semiannual date
halfyear_end <- function(date) {
  day_lookup <- c(30,30,30,30,30,30,31,31,31,31,31,31)
  months <- (floor((month(date)-1)/6)+1)*6
  days <- day_lookup[months]
  as.Date(ISOdate(year(date),months,days))
}

# Return end of quarter date
quarter_end <- function(date) {
  day_lookup <- c(31,31,31,30,30,30,30,30,30,31,31,31)
  months <- (floor((month(date)-1)/3)+1)*3
  days <- day_lookup[months]
  as.Date(ISOdate(year(date),months,days))
}

# Return end of month date
month_end <- function(date) {
  day(date) <- days_in_month(date)
  date
}