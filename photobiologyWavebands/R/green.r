#' Definition of green waveband
#'
#' Green radiation according to ISO or as commonly defined in plant
#' photobiology, no weighting applied.
#'
#' @param std a character string "ISO", "Sellaro" or "Plant"
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @note When released, this package will replace the package UVcalc.
#' @examples
#' Green()
#' Green("ISO") # 500 to 570
#' Green("Sellaro") # 500 to 570 nm
#'
#' @family unweighted wavebands
#'
Green <- function(std="ISO"){
  label="Green"
  if (std=="ISO") {
    return(new_waveband(500, 570, wb.name=paste("Green", std, sep="."), wb.label=label))
  } else if (std=="Sellaro"){
    return(new_waveband(500, 570, wb.name=paste("Green", std, sep="."), wb.label=label))
  } else {
    warning("'std' argument value not implemented.")
    return(NA)
  }
}
