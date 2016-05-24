#' Converts date from standard notation to week notation according to ISO 8601
#'
#' This function returns the year, the week of the year, and the day of week of a given date according to ISO 8601.
#' It is an substitute for the \code{\%Y-W\%V-\%u} format which is not implemented on Windows.
#'
#' According to ISO 8601, the year of the week can differ from the calendar year (see the examples).
#'
#' @param date Vector which can be coerced to class \code{Date}
#' @return A character vector of year, week, and weekday in format "\code{\%Y-W\%V-\%u}"
#' @seealso \code{\link{strptime}} for a description of the date formats and references on ISO 8601. 
#' @export
#' @author Uwe Block \email{u.block.mz@@googlemail.com}
#' @examples
#' x <- paste(1999:2011, "-12-31", sep = "")
#' y <- as.Date(x)
#' data.frame(date = format(y), weekdate = date2ISOweek(y))
#' data.frame(date = x, weekdate = date2ISOweek(x))
date2ISOweek <- function(date) {
  return(ifelse(is.na(date), NA_character_, paste(ISOweek(date), ISOweekday(date), sep = "-")))
}