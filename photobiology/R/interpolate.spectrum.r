#' Calculate spectral values at a different set of wavelengths
#'
#' For example interpolate spectral irradiance (or spectral transmittance)
#' values at new wavelengths values.
#'
#' @param w.length.in numeric array of wavelengths (nm)
#' @param s.irrad a numeric array of spectral values
#' @param w.length.out numeric array of wavelengths (nm)
#' @param fill a value to be assigned to out of range wavelengths
#'
#' @return a numeric array of interpolated spectral values
#'
#' @export
#' @note The current version of interpolate uses \code{spline} if fewer than 25
#' data points are available. Otherwise it uses \code{approx}. In the first case
#' a cubic spline is used, in the second case linear interpolation, which should
#' be faster.
#'
#' @examples
#' 
#' my.w.length <- 300:700
#' my.s.e.irrad <-
#'   with(sun.data, interpolate_spectrum(w.length, s.e.irrad, my.w.length))
#' plot(my.s.e.irrad ~ my.w.length)
#' lines(s.e.irrad ~ w.length, data=sun.data)
#'
interpolate_spectrum <- function(w.length.in, s.irrad, w.length.out, fill=NA) {
  if (is.null(fill) && (w.length.out[1] < w.length.in[1] ||
                                w.length.out[length(w.length.out)] > w.length.in[length(w.length.in)])) {
    stop("Extrapolation attempted with fill=NULL")
  }
  selector <- w.length.out >= w.length.in[1] & w.length.out <= w.length.in[length(w.length.in)]
  s.irrad.out <- numeric(length(w.length.out))
  if (!is.null(fill)){
    s.irrad.out[!selector] <- fill
  }
  if (sum(selector) < 1) {
    NULL
  } else if (sum(selector) <= 25) {
    s.irrad.out[selector] <- stats::spline(w.length.in, s.irrad, xout=w.length.out[selector])$y
  } else {
    s.irrad.out[selector] <- stats::approx(w.length.in, s.irrad, xout=w.length.out[selector])$y
  }
  return(s.irrad.out)
}
