#' Gives values for the CH4 production from pectin BSWF as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances. The
#' returned values are on quantum based effectiveness relative units.
#'
#' @param w.length numeric array of wavelengths (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for the BSWF normalized
#' as in the original source (300 nm) and based on energy effectiveness.
#'
#' @references
#' Bloom, A. A.; Lee-Taylor, J.; Madronich, S.; Messenger, D. J.; Palmer, P. I.;
#' Reay, D. S. & McLeod, A. R. (2010) Global methane emission estimates from
#' ultraviolet irradiation of terrestrial plant foliage. New Phytologist,
#' Blackwell Publishing Ltd, 187, 417--425 .
#'
#' @export
#' @examples
#' CH4_e_fun(293:400)
#'
#' @family BSWF functions
#'
CH4_e_fun <-
function(w.length){
  10^((300 - w.length) / 80)
}

