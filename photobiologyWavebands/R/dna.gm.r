#' Definition of DNA damage (SETLOW) weighted waveband
#'
#' Naked DNA damage BSWF, Green and Miller's formulation.
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
#' @seealso \code{\link{new_waveband}}  \code{\link{waveband}}
#'
#' @examples
#' DNA_GM()
#' DNA_GM(300)
#'
#' @family BSWF weighted wavebands
#'
DNA_GM <- function(norm=300, w.low=275, w.high=400) {
  new_waveband(w.low=w.low, w.high=w.high,
               weight="SWF", SWF.q.fun=DNA_GM_q_fun, SWF.norm=300, norm=norm,
               wb.name=paste("DNA.GM", as.character(norm), sep="."), wb.label="DNA G&M")
}

#' Definition of DNA damage (SETLOW) weighted waveband
#'
#' Naked DNA damage BSWF, Green and Miller's formulation.
#'
#' @seealso \code{\link{DNA_GM}}
#'
#' @export
#'
#' @keywords internal
#'
DNA.GM <- DNA_GM
