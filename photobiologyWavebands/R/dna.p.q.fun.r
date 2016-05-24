#' Gives values for plant DNA BSWF (Quaite) as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances. It uses
#' the formulation proposed by Musil.
#'
#' @param w.length numeric array of w.length (nm)
#'
#' @return a numeric array of the same length as \code{w.length} with values for
#'   the BSWF normalized as in the original source.  The returned values are
#'   based on quantum effectiveness units.
#'
#'
#' @export
#' @examples
#' DNA_P_q_fun(293:400)
#'
#'
#' @family BSWF functions
#'
DNA_P_q_fun <-
function(w.length){
    QUAITE_MUSIL.quantum290 <- numeric(length(w.length))
    QUAITE_MUSIL.quantum290[w.length<400] <-
      22.657e-3*7.98e-16*exp(1.118e4/w.length[w.length<400])
    QUAITE_MUSIL.quantum290[w.length >= 400] <- 0
    return(QUAITE_MUSIL.quantum290)
}

