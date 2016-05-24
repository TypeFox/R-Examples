#' Gives irradiance from spectral irradiance.
#'
#' This function gives the result of integrating spectral irradiance over
#' wavelengths.
#'
#' @param x numeric array
#' @param y numeric array
#'
#' @return a single numeric value with no change in scale factor: e.g. [W m-2
#'   nm-1] -> [W m-2]
#' @export
#'
#' @examples
#' with(sun.data, integrate_xy(w.length, s.e.irrad))
#'
integrate_xy <- function(x, y) {
  j <- length(x)
  sum((y[1:j-1] + y[2:j]) * diff(x)) / 2
}
