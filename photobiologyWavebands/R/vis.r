#' Definition of VIS waveband
#'
#' Visible (to humnas) radiation (380...760 nm), no weighting
#' applied.
#'
#' @param std a character string "ISO"
#'
#' @return A waveband object wavelength defining a wavelength range.
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' VIS()
#' VIS("ISO")
#'
#' @family unweighted wavebands
#'
VIS <- function(std="ISO"){
  if (std=="ISO") {
    return(new_waveband(380, 760, wb.name="VIS.ISO", wb.label="VIS"))
  } else {
    warning("'std' argument value not implemented.")
    return(NA)
  }
}
