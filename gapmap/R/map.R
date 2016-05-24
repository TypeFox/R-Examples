#'A function to map values in a range
#'
#'This function maps a value in one range to another range.
#'
#'
#' @param value input value
#' @param start1 lower bound of the value's current range
#' @param stop1 upper bound of the value's current range
#' @param start2 lower bound of the value's taget range
#' @param stop2 upper bound of the value's target range
#' @export map
#' @aliases map
#' @return a numeric value
#' @keywords internal
#' 


map <- function(value, start1, stop1, start2, stop2) {
  return (start2 + (stop2 - start2) * ((value - start1) / (stop1 - start1)))
}


# map(1, 0, 10, 0, 10)
# map(5, 0, 10, 0, 10)
# map(10, 0, 10, 0, 10)