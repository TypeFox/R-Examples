#' Gives values for GPAS BSWF (Micheletti's formulation) as a function of
#' wavelength
#'
#' This function gives a set of numeric multipliers that can be used as a weight
#' to calculate effective doses and irradiances. The BSWF is normalized at 300
#' nm.
#'
#' @param w.length numeric array of w.length (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for
#'   the BSWF normalized as in the original source.  The returned values are
#'   based on quantum effectiveness units.
#'
#' @references [1]Caldwell, M. M. (1971) Solar UV irradiation and the growth and
#' development of higher plants. In Giese, A. C. (Ed.) Photophysiology, Academic
#' Press, 1971, 6, 131-177
#'
#' [2] Micheletti, M. I. and R. D. Piacentini (2002) Irradiancia espetral solar
#' UV-B y su relación con la efectividad de daño biológico a las plantas. ANALES
#' AFA, 13, 242-248
#'
#' [3] Micheletti, M. I.; Piacentini, R. D. & Madronich, S. (2003) Sensitivity
#' of Biologically Active UV Radiation to Stratospheric Ozone Changes: Effects
#' of Action Spectrum Shape and Wavelength Range Photochemistry and
#' Photobiology, 78, 456-461
#'
#' @note In the original publication [2] describing the formulation, the
#'   long-end wavelength boundary is not specified, but 313.3 nm is usually
#'   used. The equation is coded here with the limit at 342 nm as at longer
#'   wavelengths the values increase with increasing wavelength. The effect on
#'   the RAF and doses of changing this boundary ican be substantial, and has
#'   been analysed by Micheletti et al. [3].
#'
#'
#' @export
#' @examples
#' GEN_M_q_fun(293:400)
#'
#' @family BSWF functions
#'
GEN_M_q_fun <-
  function(w.length){
    spectral_weights <-
      570.25 - 4.70144 * w.length + 0.01274 * w.length^2 - 1.13118E-5 * w.length^3

    ifelse(spectral_weights < 0 | w.length > 342, 0, spectral_weights)

  }


