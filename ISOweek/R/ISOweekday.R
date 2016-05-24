#' Day of week according to ISO 8601
#'
#' This function returns the weekday of a given date according to ISO 8601.
#' It is an substitute for the "\code{\%u}" format which is not implemented on Windows.
#' 
#' @param date Vector which can be coerced to class \code{Date}
#' @return An integer vector of weekdays (1-7, Monday is 1)
#' @author Uwe Block \email{u.block.mz@@googlemail.com}
#' @seealso \code{\link{strptime}}
#' @export
#' @examples
#' x <- paste(1999:2011, "-12-31", sep = "")
#' y <- as.Date(x)
#' data.frame(date = format(y), weekday = ISOweekday(y))
#' data.frame(date = x, weekday = ISOweekday(x))
ISOweekday <- function(date) {
  date <- as.Date(date)
  return(as.integer((as.integer(format(date, "%w"))+6) %% 7 + 1))
}

