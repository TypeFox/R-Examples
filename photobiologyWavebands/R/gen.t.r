#' Definition of GPAS (Timijan) weighted waveband
#'
#' Generalized Plant Action BSWF of Caldwell [1] as formulated by Timijan et al.
#' [2]
#'
#' @param norm normalization wavelength (nm)
#' @param w.low short-end boundary wavelength (nm)
#' @param w.high long-end boundary wavelength (nm)
#'
#' @return a waveband object wavelength defining wavelength range, weighting
#'   function and normalization wavelength.
#'
#' @references [1] Caldwell, M. M. (1971) Solar UV irradiation and the growth
#' and development of higher plants. In Giese, A. C. (Ed.) Photophysiology,
#' Academic Press, 1971, 6, 131-177
#'
#' [2] Thimijan RW, Cams HR, Campbell L. (1978) Radiation sources and related
#' environmental control for biological and climatic eflFects of UV research.
#' Final report EPA-IAG-D6-0168. Washington: Environmental Protection Agency.
#'
#' @export
#' @seealso \code{\link{GEN.G}} \code{\link{GEN.M}} \code{\link{PG}}
#'   and \code{\link[photobiology]{waveband}}
#' @examples
#' GEN_T()
#' GEN_T(300)
#'
#' @family BSWF weighted wavebands
#'
GEN_T <- function(norm=300, w.low=275, w.high=345) {
  new_waveband(w.low=w.low, w.high=w.high, weight="SWF", SWF.q.fun=GEN_T_q_fun, SWF.norm=300,
               norm=norm, wb.name=paste("GEN.T", as.character(norm), sep="."), wb.label="GEN(T)")
}

#' Definition of GPAS (Timijan) weighted waveband
#'
#' Generalized Plant Action BSWF of Caldwell as formulated by Timijan et al.
#'
#' @export
#'
#' @seealso \code{\link{GEN_T}}
#'
#' @keywords internal
#'
GEN.T <- GEN_T
