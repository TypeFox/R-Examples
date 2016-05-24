#' Calculate deuterim lamp output spectrum from fitted constants
#'
#' @description Calculate values by means of a nth degree polynomial from
#' user-supplied constants (for example from a lamp calibartion certificate).
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
#' @export D2_spectrum
#'
#' @note This is function is valid for wavelengths in the range 180 nm to 495
#'   nm, for wavelengths outside this range NAs are returned.
#' @examples
#' D2_spectrum(200)
#' D2_spectrum(170:220)
#'
D2_spectrum <- function(w.length, k=D2.UV653, fill=NULL) {
  pws <- (length(k) - 1):0
  fill.selector <- w.length < 190 | w.length > 450
  if (is.null(fill)) {
    w.length <- w.length[!fill.selector]
    fill.selector <- rep(FALSE, length(w.length))
    indexes <- 1:length(w.length)
  } else {
    indexes <- which(w.length >= 190 & w.length <= 450)
  }
  s.e.irrad <- numeric(length(w.length))
  for (i in indexes) {
    s.e.irrad[i] <- sum(w.length[i]^pws * k)
  }
  s.e.irrad[fill.selector] <- fill
  out.data <- source_spct(w.length, s.e.irrad)
  comment(out.data) <- paste("Fitted spectrum for:", comment(k))
  return(out.data)
}

D2.UV653 <- c(-4.1090384E-17, 7.6305376E-14, -5.6315241E-11, 2.0714647E-08, -3.8162638E-06, 0.0002842542)
comment(D2.UV653) <- "D2 lamp UV-653"

D2.UV586 <- c(7.1215397E-18, -4.0918226E-15, -2.9069045E-12, 3.0340054E-09, -8.9554589E-07, 9.1200585E-05)
comment(D2.UV586) <- "D2 lamp UV-586"

D2.UV654 <- c(-6.2181585E-17, 1.135044E-13, -8.2196585E-11, 2.9598246E-08, -5.3217495E-06, 0.00038528517)
comment(D2.UV654) <- "D2 lamp UV-654"



