#' Determine numeric age
#' 
#' Determines the numeric age in years from a date range.
#' 
#' @param from character, beginning of the date range in \code{YYYY-MM-DD} format.
#' @param to character, end of the date range in \code{YYYY-MM-DD} format.
#' @details Helper for LexisPlotR. The numeric age gets rounded to 5 digits.
#' @return Numeric age in years.
#' @author Philipp Ottolinger
#' @export how.old
#' @examples
#' library(LexisPlotR)
#' how.old(from = "1900-01-01", to = "1905-01-01")

how.old <- function(from, to) {
  from <- as.Date(from)
  to <- as.Date(to)
  age <- round(as.numeric(to-from)/365.25, digits = 5)
  return(age)
}