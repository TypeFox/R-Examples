#' Definition of CIE weighted waveband
#'
#' Erythema BSWF
#'
#' @param norm normalization wavelength (nm)
#' @param w.low short-end boundary wavelength (nm)
#' @param w.high long-end boundary wavelength (nm)
#'
#' @return a waveband object wavelength defining wavelength range, weighting
#'   function and normalization wavelength.
#'
#' @references
#' Webb, A. (20XX)
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' CIE()
#' CIE(300)
#'
#' @family BSWF weighted wavebands
#'
CIE <- function(norm = 298, w.low = 250, w.high = 400) {
  new_waveband(w.low = w.low, w.high = w.high,
               weight = "SWF", SWF.e.fun = CIE_e_fun, SWF.norm = 298,
               norm = norm, hinges=c(250 - 1e-12, 250, 298, 328, 400 - 1e-12, 400),
               wb.name = paste("CIE98", as.character(norm), sep = "."), wb.label = "CIE98")
}
