#' Create thermal tolerance curves.
#'
#' Creates a quadratic thermal tolerance curve of the form: width_param * (temp
#' - optim_temp)^2 + max_a Negative values are *not* returned as 0 for speed of
#' computation. You should check for this after.
#'
#' @param temp The input temperature value.
#' @param optim_temp The optimal temperature.
#' @param max_a The maximum productivity parameter `a` from a Ricker model (or
#'   whatever the y-axis value is you want to return).
#' @param width_param A parameter to control the width of the parabola. Smaller
#'   numbers make wider parabolas.
#' @return A productivity parameter given the location on a thermal tolerance
#'   curve.
#' @export
#' @examples
#' x <- seq(5, 30, length.out = 200)
#' plot(x, thermal_curve_a(x), ylab = "a", xlab = "Temperature", type
#' = "l")

thermal_curve_a <- function(temp, optim_temp = 15, max_a = 1.4, width_param = 0.02) {
  a <- -width_param * (temp - optim_temp)^2 + max_a
  #return(ifelse(a < 0, 0, a))
  a
}

