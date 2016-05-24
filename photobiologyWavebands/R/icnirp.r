#' Definition of ICNIRP 2004 weighted waveband
#'
#' ICNIRP 2004 BSWF waveband constructor. This BSWF is used for the
#' determination of exposure limits (EL) for workers, and includes a safety margin
#' as it is based on eye and the non-pathologic response of the most sensitive
#' human skin types when not tanned.
#'
#' @param norm normalization wavelength (nm)
#' @param w.low short-end boundary wavelength (nm)
#' @param w.high long-end boundary wavelength (nm)
#'
#' @return a waveband object defining wavelength range, weighting function
#' and normalization wavelength.
#'
#' @references
#'   INTERNATIONAL COMMISSION ON NON-IONIZING RADIATION PROTECTION (2004) ICNIRP
#'   GUIDELINES ON LIMITS OF EXPOSURE TO ULTRAVIOLET RADIATION OF WAVELENGTHS
#'   BETWEEN 180 nm AND 400 nm (INCOHERENT OPTICAL RADIATION). HEALTH PHYSICS
#'   87(2):171-186.
#'   \url{http://www.icnirp.org/cms/upload/publications/ICNIRPUV2004.pdf}
#'
#'
#' @export
#'
#' @seealso \code{\link{new_waveband}}  \code{\link{waveband}}
#'
#' @examples
#' ICNIRP()
#'
#' @family BSWF weighted wavebands
#'
ICNIRP <- function(norm = 270, w.low = 210, w.high = 400) {
  new_waveband(w.low = w.low, w.high = w.high,
               weight = "SWF", SWF.e.fun = ICNIRP_e_fun, SWF.norm = 270,
               norm = norm, hinges = c(210 - 1e-12, 210, 270, 270 + 1e-12, 300, 300 + 1e-12, 400, 400 + 1e-12),
               wb.name = paste("ICNIRP", as.character(norm), sep = "."), wb.label = "ICNIRP")
}
