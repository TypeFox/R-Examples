#' @title Four parameter logistic desirability function
#'
#' @description Maps a numeric variable to a 0-1 scale with a logistic function.
#'
#' @details This function uses a four parameter logistic model to map
#' a numeric variable onto a 0-1 scale. Whether high or low values are
#' deemed desirable can be controlled with the \code{hill} parameter;
#' when \code{hill} > 0 high values are desirable and when \code{hill}
#' < 0 low values are desirable
#'
#' Note that if the data contain both positive and negative values
#' this function does not provide a monotonic mapping (see example).
#'
#' @param x Vector of numeric or integer values.
#' @param des.min,des.max The lower and upper asymptotes of the
#' function. Defaults to zero and one, respectively.
#' @param hill Hill coefficient. It controls the steepness and direction
#' of the slope. A value greater than zero has a positive slope and
#' a value less than zero has a negative slope. The higher the absolute
#' value, the steeper the slope.
#' @param inflec Inflection point. Is the point on the x-axis where
#' the curvature of the function changes from concave upwards to
#' concave downwards (or vice versa).
#'
#' @return Numeric vector of desirability values.
#' @seealso \code{\link{d.low}},  \code{\link{d.high}}
#' 
#' @examples
#' # High values are desirable
#' x1 <- seq(80, 120, 0.01)
#' d1 <- d.4pl(x = x1, hill = 20, inflec = 100)
#' plot(d1 ~ x1, type="l")
#'
#' # Low values are desirable (negative slope), with a minimum
#' # desirability of 0.3
#' d2 <- d.4pl(x = x1, hill = -30, inflec = 100, des.min=0.3)
#' plot(d2 ~ x1, type="l", ylim=c(0,1))
#'
#' # Beware of how the function behaves when the data contain both
#' # positive and negative values
#' x2 <- seq(-20, 20, 0.01)
#' d3 <- d.4pl(x = x2, hill = 20, inflec = 1)
#' plot(d3 ~ x2, type="l")

d.4pl <- function(x, hill, inflec, des.min = 0, des.max = 1){

  if (hill == 0) stop("The Hill coefficient must not equal zero\n")
  if (des.min < 0 | des.min > 1) stop("des.min must be between zero and one\n")
  if (des.max < 0 | des.max > 1) stop("des.max must be between zero and one\n")
  
  y <- ((des.min - des.max) / (1 + ((x/inflec)^hill))) + des.max
  
  return(y)
}
