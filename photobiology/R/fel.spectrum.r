#' Incandescent "FEL" lamp emission spectrum
#'
#' @description Calculate values by means of a nth degree polynomial from
#' user-supplied constants (for example from a lamp calibration certificate).
#'
#' @param w.length numeric vector of wavelengths (nm) for output
#' @param k a numeric vector with n constants for the polynomial
#' @param fill if NA, no extrapolation is done, and NA is returned for
#'   wavelengths outside the range of the input. If NULL then the tails are
#'   deleted. If 0 then the tails are set to zero, etc. NULL is default.
#'
#' @return a dataframe with four numeric vectors with wavelength values
#'   (w.length), energy and photon irradiance (s.e.irrad, s.q.irrad) depending
#'   on the argument passed to unit.out (s.irrad).
#'
#' @export FEL_spectrum
#'
#' @note This is function is valid for wavelengths in the range 180 nm to 495
#'   nm, for wavelengths outside this range NAs are returned.
#' @examples
#' FEL_spectrum(200)
#' FEL_spectrum(170:220)

FEL_spectrum <- function(w.length, k=FEL, fill=NA) {
  pws <- (length(k$kb) - 1):0
  fill.selector <- w.length < 250 | w.length > 900
  if (is.null(fill)) {
    w.length <- w.length[!fill.selector]
    fill.selector <- rep(FALSE, length(w.length))
    indexes <- 1:length(w.length)
  } else {
    indexes <- which(w.length >= 250 & w.length <= 900)
  }
  s.e.irrad <- numeric(length(w.length))
  for (i in indexes) {
    s.e.irrad[i] <- sum((k$kb * w.length[i]^pws) * k$kc /
                          ((w.length[i] * 1e-9)^5 * (exp(0.014388 / (w.length[i] * 1e-9) / k$TK) - 1)))
  }
  s.e.irrad[fill.selector] <- fill
  out.data <- source_spct(w.length, s.e.irrad)
  comment(out.data) <- paste("Fitted spectrum for:", comment(k))
  return(out.data)
}

FEL <- list(TK = 3106.87561155307, kc = 1.5032470318064E-10 * 3.74183488349896E-16 * 0.0000001,
            kb = c(1.48299373059776E-19, -6.73783136663075E-16, 1.28071973512943E-12, -1.31864190324975E-09, 7.93362535586965E-07, -0.000278602532732717, 0.052698171960183, -3.66729840729808))
comment(FEL) <- "FEL lamp BN-9101-165, from calibration document T-R_387.xls"

