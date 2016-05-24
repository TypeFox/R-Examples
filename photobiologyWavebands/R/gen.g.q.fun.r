#' Gives values for GPAS BSWF (Green's formulation) as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances. The BSWF is normalized
#' at 280 nm.
#'
#' @param w.length numeric array of w.length (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for the BSWF normalized
#' as in the original source.  The returned values are based on quantum effectiveness units.
#'
#' @references
#' [1] Caldwell, M. M. (1971) Solar UV irradiation and the growth and development
#' of higher plants. In Giese, A. C. (Ed.) Photophysiology, Academic Press,
#' 1971, 6, 131-177
#'
#' [2] Green, A. E. S.; Sawada, T. & Shettle, E. P. (1974) The middle
#' ultraviolet reaching the ground Photochemistry and Photobiology, 1974, 19,
#' 251-259
#'
#' [3] Micheletti, M. I.; Piacentini, R. D. & Madronich, S. (2003) Sensitivity
#' of Biologically Active UV Radiation to Stratospheric Ozone Changes: Effects
#' of Action Spectrum Shape
#' and Wavelength Range Photochemistry and Photobiology, 78, 456-461
#'
#' @note In the original publication [2] describing the formulation, the long-end
#' wavelength boundary is specified as 313.3 nm. The equation
#' is coded here with no such limit so that any limit can be set when defining the
#' waveband. We do so because in some cases it is of interest to vary this limit
#' in sensitivity analyses.
#' The effect on the RAF and doses of changing this boundary is substantial, and
#' has been analysed by Micheletti et al. [3].
#'
#'
#' @export
#' @examples
#' GEN_G_q_fun(293:400)
#'
#' @family BSWF functions
#'
GEN_G_q_fun <-
  function(w.length){
    spectral_weights <-
      2.618 * (1.0 - (w.length / 313.3)^2) * exp((300 - w.length) / 31.08)

    ifelse(spectral_weights < 0, 0, spectral_weights)
  }
