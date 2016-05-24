#' Gives values for naked DNA BSWF (SETLOW) as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used as a weight
#' to calculate effective doses and irradiances. It uses the seldom used Green
#' and Miller formulation.
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
#' DNA_GM_q_fun(293:400)
#'
#' @family BSWF functions
#'
DNA_GM_q_fun <-
function(w.length){
    SETLOW_GM.quantum300 <- numeric(length(w.length))
    SETLOW_GM.quantum300 <- exp(13.82*(1/(1+exp((w.length-310)/9))-1))*30.675
    return(SETLOW_GM.quantum300)
}

