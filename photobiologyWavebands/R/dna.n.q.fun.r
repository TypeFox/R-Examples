#' Gives values for naked DNA BSWF (SETLOW) as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances.
#'
#' @param w.length numeric array of w.length (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for
#'   the BSWF normalized as in the original source.  The returned values are
#'   based on quantum effectiveness units.
#'
#' @note The digitized data as used in the TUV model covers the wavelength range
#'   from 256 nm to 364 nm. For longer wavelengths we set the value to zero, and
#'   for shorter wavelengths we extrapolate the value for 256 nm.
#'
#'
#' @export
#' @examples
#' DNA_N_q_fun(293:400)
DNA_N_q_fun <-
  function(w.length){
    wl.within <- w.length >= 256 & w.length <= 364
    spectral_weights <- numeric(length(w.length))
    spectral_weights[w.length < 256] <- NA # the value at 256 nm
    if (any(wl.within)) { # avoids error in spline when xout is empty
      spectral_weights[wl.within] <-
        stats::spline(photobiologyWavebands::SetlowTUV.spct$w.length,
                      photobiologyWavebands::SetlowTUV.spct$s.q.response,
                      xout = w.length[wl.within])$y
    }
    spectral_weights[w.length > 364] <- 0.0
    return(spectral_weights)
  }

