#' Definition of FR waveband
#'
#' Far-red radiation according to "ISO" (not defined) or as commonly defined in
#' plant photobiology, "Smith10" (725-735 nm), "Smith20" (720-740 nm), "Inada"
#' (700-800 nm), "Warrington" (700-850 nm), and "Sellaro" (700-750 nm), and
#' "BTV" (700-770 nm), as defined in a recent handbook. No weighting applied.
#'
#' @param std a character string, defaults to "ISO", as for other colour
#'   definitions, which in this case returns \code{NA}.
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' Far_red()
#' Far_red("ISO")
#' Far_red("Smith")
#' Far_red("BTV")
#'
#' @family unweighted wavebands
#'
Far_red <- function(std="ISO"){
  label="FR"
  if (std=="Smith") {
    warning("The definition of 'Smith' defaults to 'Smith10', to restore old behaviour use 'Smith20'.")
    std <- "Smith10"
  }
  if (std=="ISO") {
    warning("'ISO gives no standard definition of far-red.")
    return(NA)
  } else if (std=="Smith20"){
    return(new_waveband(720, 740, wb.name=paste("FarRed", std, sep="."), wb.label=label))
  } else if (std=="Smith10"){
    return(new_waveband(725, 735, wb.name=paste("FarRed", std, sep="."), wb.label=label))
  } else if (std=="Inada"){
    return(new_waveband(700, 800, wb.name=paste("FarRed", std, sep="."), wb.label=label))
  } else if (std=="Warrington"){
    return(new_waveband(700, 850, wb.name=paste("FarRed", std, sep="."), wb.label=label))
  } else if (std=="Sellaro"){
    return(new_waveband(700, 750, wb.name=paste("FarRed", std, sep="."), wb.label=label))
  } else {
    warning("'std' argument value not implemented.")
    return(NA)
  }
}
