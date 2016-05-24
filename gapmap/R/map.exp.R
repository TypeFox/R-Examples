#'A function to map values in a range exponentially
#'
#'This function maps a value in one range to another range exponentially.
#'
#'
#' @param value input value
#' @param start1 lower bound of the value's current range
#' @param stop1 upper bound of the value's current range
#' @param start2 lower bound of the value's taget range
#' @param stop2 upper bound of the value's target range
#' @param scale scale log base
#' @export map.exp
#' @aliases map.exp
#' @return a numeric value
#' @keywords internal
#' 


map.exp <- function(value, start1, stop1, start2, stop2, scale = 0.5) {
  f <- (value - start1) / ( stop1 - start1)
  flog <- f ^ (1/scale)
  return (flog*(stop2- start2))
}

# map.exp(1, 0, 10, 0, 1, base= 0.5)
# map.exp(5, 0, 10, 0, 1, base= 0.5)
# map.exp(10, 0, 10, 0, 1, base= 0.5)



