#' Gives values for the Plant Growth BSWF as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances. The
#' returned values are on quantum based effectiveness relative units.
#'
#' @param w.length numeric array of wavelengths (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for
#'   the BSWF normalized as in the original source (300 nm)
#'
#' @note We follow the original defition here for the equation, with no
#' limtation to the wavelength range. However, be aware that in practice
#' it is not used for long wavelengths (different limits between 366 nm and
#' 400 nm have been used by different authors).
#'
#' @export
#' @examples
#' PG_q_fun(293:400)
#'
#' @family BSWF functions
#'
PG_q_fun <-
  function(w.length){
      exp(4.688272*exp(-exp(0.1703411*(w.length-307.867)/1.15)) +
            ((390-w.length)/121.7557-4.183832))
  }

