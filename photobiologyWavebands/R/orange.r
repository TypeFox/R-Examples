#' Definition of orange waveband
#'
#' Orange radiation (591...610 nm), no weighting
#' applied.
#'
#' @param std a character string "ISO"
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @export
#' @examples
#' Orange()
#' Orange("ISO")
#'
#' @family unweighted wavebands
#'
Orange <- function(std="ISO"){
  if (std=="ISO") {
    return(new_waveband(591, 610, wb.name="Orange.ISO", wb.label="Orange"))
  } else {
    warning("'std' argument value not implemented.")
    return(NA)
  }
}
