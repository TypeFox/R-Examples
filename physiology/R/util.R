#' @title age from birth and reference dates
#' @param birth.date Date of birth, either as a \code{Date} or something which
#'   will be converted to a \code{Date}
#' @param ref.date Date at which to calculate age, defaults to current date,
#'   either as a \code{Date} or something which will be converted to a
#'   \code{Date}
#' @param unit character of length, one of "year" or "day".
#' @return integer vector
#' @examples
#' age_from_dates("2014-11-08", "2014-12-31", unit = "day")
#' age_from_dates("1981-07-09", "2014-06-29", unit = "year")
#' @export
age_from_dates <- function(birth.date, ref.date = Sys.Date(),
                         unit = c("year", "day")) {
  unit <- match.arg(unit)
  birth.date <- as.Date(birth.date)
  ref.date <- as.Date(ref.date)
  pdt <- as.POSIXlt(c(birth.date, ref.date))

  years <- pdt$year[2] - pdt$year[1]
  months <- pdt$mon[2] - pdt$mon[1] # of year
  days <- pdt$mday[2] - pdt$mday[1] # of month
  year.correct <- (months + days / 32) < 0
  if (unit == "year")
    return(years - year.correct)  # not birthday yet this year

  as.integer(ref.date - birth.date) # days
}

# TODO: print class to show days to 30d, months to 24 months, then years
# automatically
