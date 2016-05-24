#' Gives values for GPAS BSWF (Timijan's formulation) as a function of
#' wavelength
#'
#' This function gives a set of numeric multipliers that can be used as a weight
#' to calculate effective doses and irradiances.
#'
#' @param w.length numeric array of w.length (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for
#'   the BSWF normalized as in the original source.  The returned values are
#'   based on quantum effectiveness units.
#'
#' @note For wavelengths shorter than 256 nm the value returned by the equation
#'   starts decreasing, but we instead extrapolate this maximum value, obtained
#'   at 256 nm, to shorter wavelengths. For wavelengths longer than 345 nm we
#'   return zero, as is usual parctice.
#'
#'
#' @export
#' @examples
#' GEN_T_q_fun(293:400)
#'
#' @family BSWF functions
#'
GEN_T_q_fun <-
function(w.length){
    wl.within <- w.length >= 265 & w.length <= 345
    spectral_weights <- numeric(length(w.length))
    spectral_weights[w.length < 265] <- 16.08324
    if (any(wl.within)) {
      spectral_weights[wl.within] <-
      exp(-(((265-w.length[wl.within])/21)^2))/0.06217653
    }
    spectral_weights[w.length > 345] <- 0.0
    return(spectral_weights)
}

