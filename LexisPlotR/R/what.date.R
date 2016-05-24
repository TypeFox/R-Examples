#' Determine the date a certain age is reached
#' 
#' Determines the date a certain age is reached given an earlier date.
#' @param date character, set the reference date in \code{YYYY-MM-DD} format.
#' @param age numeric, set an age to be reached.
#' @details Helper for LexisPlotR.
#' @return The date \code{age} is reached when counting from \code{date}.
#' @author Philipp Ottolinger
#' @export what.date
#' @examples 
#' library(LexisPlotR)
#' what.date(date = "1900-01-01", age = 3)

what.date <- function(date, age) {
  if (!is.numeric(age)) { stop("No numeric age.") }
  date <- as.Date(date)
  date <- date + age*365.25
  return(date)
}