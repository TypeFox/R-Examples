#' Gives values for FLAV BSWF (flavonoid) as a function of wavelength
#'
#' This function gives a set of numeric multipliers that can be used
#' as a weight to calculate effective doses and irradiances. It is the
#' action spectrum for the accumulation of mesembryanthin.
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
#' FLAV_q_fun(293:400)
#'
#' @family BSWF functions
#'
FLAV_q_fun <-
function(w.length){
    FLAV.quantum <- numeric(length(w.length))
    FLAV.quantum[w.length >= 280 & w.length <= 346] <-
      exp(45.0 - 0.15 * w.length[w.length >= 280 & w.length <= 346])
    FLAV.quantum[w.length > 346] <- 0.0
    FLAV.quantum[w.length < 280] <- exp(45.0 - 0.15 * 280)
    if (w.length[1] < 280) {
      warning("FLAV BSWF is extrapolated for w.length < 280 nm\n to the value for 280 nm")
    }
    return(FLAV.quantum)
}

