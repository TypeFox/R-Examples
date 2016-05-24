#' Weekday as integer number (0-6, Monday is 0) 
#'
#' This internal function returns the weekday of a given date.
#' 
#' The week starts on Monday and ends on Sunday.
#'
#' @param date Vector which can be coerced to class \code{Date}
#' @return An integer vector of weekdays (0-6, Monday is 0)
#' @seealso \code{\link{ISOweekday}}
#' @keywords internal
weekday0 <- function(date) {
  return(ISOweekday(date) - 1L)
}

#' Date of the nearest Thursday of a given date
#'
#' This internal function returns the date of the Thursday of the week in which the given date is located.
#' 
#' The week starts on Monday and ends on Sunday.
#'
#' @param date Vector which can be coerced to class \code{Date}
#' @return A vector of dates of the nearest Thursdays
#' @keywords internal
thursday0 <- function(date) {
  date <- as.Date(date)
  return(date - weekday0(date) + 3)
}

#' Calendar year of a given date
#'
#' This internal function returns the year with century as integer.
#'
#' @param date Vector which can be coerced to class \code{Date}
#' @return An integer vector of years
#' @keywords internal
year0 <- function(date) {
  date <- as.Date(date)
  return(as.integer(format(date, "%Y")))
}
