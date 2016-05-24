#' Definition of DNA damage (SETLOW) weighted waveband
#'
#' Naked DNA damage BSWF
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
#' DNA_N()
#' DNA_N(300)
#'
#' @family BSWF weighted wavebands
#'
DNA_N <- function(norm=300, w.low=275, w.high=400) {
  new_waveband(w.low=w.low, w.high=w.high, SWF.norm=300,
               weight="SWF", SWF.q.fun=DNA_N_q_fun,
               norm=norm, wb.name=paste("DNA.N", as.character(norm), sep="."), wb.label="DNA Naked")
}

#' Definition of DNA damage (SETLOW) weighted waveband
#'
#' Naked DNA damage BSWF
#'
#' @export
#'
#' @seealso \code{\link{DNA_N}}
#'
#' @keywords internal
#'
DNA.N <- DNA_N
