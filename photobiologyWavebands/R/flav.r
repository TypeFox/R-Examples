#' Definition of FLAV BSWF flavonoids
#'
#' Mesembryanthin accumulation BSWF, data and formulation from Ibdah et al.
#'
#' @param norm normalization wavelength (nm)
#' @param w.low short-end boundary wavelength (nm)
#' @param w.high long-end boundary wavelength (nm)
#'
#' @return a waveband object wavelength defining wavelength range, weighting function
#' and normalization wavelength.
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' FLAV()
#' FLAV(300)
#'
#' @family BSWF weighted wavebands
#'
FLAV <- function(norm=300, w.low=275, w.high=346) {
  new_waveband(w.low=w.low, w.high=w.high,
               weight="SWF", SWF.q.fun=FLAV_q_fun, SWF.norm=300,
               norm=norm, hinges=c(280 - 1e-12,280,346 - 1e-12,346),
               wb.name=paste("FLAV", as.character(norm), sep="."), wb.label="FLAV")
}
