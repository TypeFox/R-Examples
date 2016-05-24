#' Definition of list of VIS wavebands
#'
#' Defined according to "ISO".
#'
#' @param std a character string "ISO", or "Sellaro"
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @export
#'
#' @seealso \code{\link{new_waveband}}  \code{\link{waveband}}
#'
#' @examples
#' Blue()
#' Blue("ISO")
#' Blue("Sellaro")
#'
#' @family unweighted wavebands
#'
Blue <- function(std="ISO"){
  label="Blue"
  if (std=="ISO") {
    return(new_waveband(450,500, wb.name=paste("Blue", std, sep="."), wb.label=label))
  } else if (std=="Sellaro"){
    return(new_waveband(420,490, wb.name=paste("Blue", std, sep="."), wb.label=label))
  } else {
    warning("'std' argument value not implemented.")
    return(NA)
  }
}
