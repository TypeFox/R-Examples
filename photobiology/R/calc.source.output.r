#' Light-source spectral output
#'
#' @description Calculate interpolated values by interpolation from
#' user-supplied spectral emission data or by name for light source data
#' included in the packages photobiologySun, photobiologyLamps, or
#' photobiologyLEDs, scaling the values.
#'
#' @param w.length.out numeric vector of wavelengths (nm) for output
#' @param w.length.in numeric vector of wavelengths (nm) for input
#' @param s.irrad.in numeric vector of spectral transmittance value (fractions
#'   or percent)
#' @param unit.in a character string "energy" or "photon"
#' @param scaled NULL, "peak", "area"; div ignored if !is.null(scaled)
#' @param fill if NA, no extrapolation is done, and NA is returned for
#'   wavelengths outside the range of the input. If NULL then the tails are
#'   deleted. If 0 then the tails are set to zero.
#'
#' @return a source_spct with three numeric vectors with wavelength values
#'   (w.length), scaled and interpolated spectral energy irradiance (s.e.irrad),
#'   scaled and interpolated spectral photon irradiance values (s.q.irrad).
#'
#' @export
#'
#' @note This is a convenience function that adds no new functionality but makes
#'   it a little easier to plot lamp spectral emission data consistently. It
#'   automates interpolation, extrapolation/trimming and scaling.
#' @examples
#' with(sun.data, calc_source_output(290:1100, w.length.in=w.length, s.irrad.in=s.e.irrad))
#'
calc_source_output <- function(w.length.out,
                               w.length.in, s.irrad.in,
                               unit.in="energy",
                               scaled=NULL, fill=NA) {

  if (!check_spectrum(w.length.in, s.irrad.in)) {
      return(NA)
  }

  # we interpolate using a spline or linear interpolation
  out.fill.selector <- w.length.out < w.length.in[1] | w.length.out > w.length.in[length(w.length.in)]
  if (is.null(fill)) {
    w.length.out <- w.length.out[!out.fill.selector]
    out.fill.selector <- rep(FALSE, length(w.length.out))
  }
  s.irrad.out <- numeric(length(w.length.out))

  if (length(w.length.out) < 25) {
    # cubic spline
    s.irrad.out[!out.fill.selector] <-
      stats::spline(w.length.in, s.irrad.in, xout=w.length.out[!out.fill.selector])$y
  } else {
    # linear interpolation
    s.irrad.out[!out.fill.selector] <-
      stats::approx(x = w.length.in, y = s.irrad.in,
             xout = w.length.out[!out.fill.selector], ties = "ordered")$y
  }

  # we check unit.in and and convert the output spectrum accordingly

  if (unit.in == "energy") {
    out.data <- e2q(source_spct(w.length = w.length.out,
                            s.e.irrad = s.irrad.out),
                    action = "add")
  } else if (unit.in == "photon") {
    out.data <- q2e(source_spct(w.length = w.length.out,
                            s.q.irrad = s.irrad.out),
                    action = "add")
  } else {
    warning("Bad argument for unit.in: ", unit.in)
    return(NA)
  }

  # do scaling

  if (!is.null(scaled)) {
    if (scaled == "peak") {
      e.div <- max(out.data$s.e.irrad, na.rm=TRUE)
      q.div <- max(out.data$s.q.irrad, na.rm=TRUE)
    } else if (scaled == "area") {
      s.irrad.na.sub <- out.data$s.e.irrad
      s.irrad.na.sub[is.na(s.irrad.na.sub)] <- 0.0
      e.div <- integrate_xy(w.length.out, s.irrad.na.sub)
      s.irrad.na.sub <- out.data$s.q.irrad
      s.irrad.na.sub[is.na(s.irrad.na.sub)] <- 0.0
      q.div <- integrate_xy(w.length.out, s.irrad.na.sub)
    } else {
      warning("Ignoring unsupported scaled argument: ", scaled)
      e.div <- q.div <- 1.0
    }
    out.data[!out.fill.selector, "s.e.irrad"] <- out.data[!out.fill.selector, "s.e.irrad"] / e.div
    out.data[!out.fill.selector, "s.q.irrad"] <- out.data[!out.fill.selector, "s.q.irrad"] / q.div
  }
  out.data[out.fill.selector, "s.e.irrad"] <- fill
  out.data[out.fill.selector, "s.q.irrad"] <- fill

  return(out.data)
}
