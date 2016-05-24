#' Definition of GPAS (Green) weighted waveband
#'
#' Generalized Plant Action BSWF of Caldwell as formulated by Green et al.
#'
#' @param norm normalization wavelength (nm)
#' @param w.low short-end boundary wavelength (nm)
#' @param w.high long-end boundary wavelength (nm)
#'
#' @return a waveband object wavelength defining wavelength range, weighting function
#' and normalization wavelength.
#'
#' @references
#' [1]Caldwell, M. M. (1971) Solar UV irradiation and the growth and development
#' of higher plants. In Giese, A. C. (Ed.) Photophysiology, Academic Press,
#' 1971, 6, 131-177
#'
#' [2] Green, A. E. S.; Sawada, T. & Shettle, E. P. (1974) The middle
#' ultraviolet reaching the ground Photochemistry and Photobiology, 1974, 19,
#' 251-259
#'
#' [3] Micheletti, M. I.; Piacentini, R. D. & Madronich, S. (2003) Sensitivity
#' of Biologically Active UV Radiation to Stratospheric Ozone Changes: Effects
#' of Action Spectrum Shape and Wavelength Range Photochemistry and
#' Photobiology, 78, 456-461
#'
#' @note In the original publication [2] describing the formulation, the long-end
#' wavelength boundary is specified as 313.3 nm. This is the default used here.
#' However, in some cases it is of interest to vary this limit in sensitivity analyses.
#' The effect on the RAF and doses of changing this boundary is substantial, and
#' has been analysed by Micheletti et al. [3].
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' GEN_G()
#' GEN_G(300)
#'
#' @family BSWF weighted wavebands
#'
GEN_G <- function(norm=300, w.low=275, w.high=313.3) {
  new_waveband(w.low=w.low, w.high=w.high, weight="SWF", SWF.q.fun=GEN_G_q_fun, SWF.norm=280,
               norm=norm, wb.name=paste("GEN.G", as.character(norm), sep="."), wb.label="GEN(G)")
}

#' Definition of GPAS (Green) weighted waveband
#'
#' Generalized Plant Action BSWF of Caldwell as formulated by Green et al.
#'
#' @export
#'
#' @seealso \code{\link{GEN_G}}
#'
#' @keywords internal
#'
GEN.G <- GEN_G
