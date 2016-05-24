#' Gives values for the erythemal BSWF as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used as a weight
#' to calculate effective doses and irradiances. The returned values are on
#' quantum based effectiveness relative units.
#'
#' @param w.length numeric array of wavelengths (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for
#'   the BSWF normalized as in the original source (298 nm) and based on quantum
#'   effectiveness.
#'
#'
#' @export
#' @examples
#' CIE_q_fun(293:400)
#'
#'
#' @family BSWF functions
#'
CIE_q_fun <-
function(w.length){
    CIE_e_fun(w.length) * 298.0 / w.length
}

