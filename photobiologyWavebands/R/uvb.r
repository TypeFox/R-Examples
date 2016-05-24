#' Definition of UV-B waveband
#'
#' UV-B according to CIE and ISO standrads: 280--315 nm. UV-B according to
#' common non-standard practice: 280--320 nm. UV-B according to medical or
#' dermatological non-standard practice: 280--320 nm.
#'
#' @param std a character string "CIE", "ISO", "medical" or "none"
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' UVB()
#' UVB("ISO")
#' UVB("CIE")
#' UVB("none")
#' UVB("medical")
#'
#' @family unweighted wavebands
#'
UVB <- function(std="ISO") {
  label <- "UVB"
  if (std=="ISO" || std=="CIE"){
    return(new_waveband(w.low=280, w.high=315,
                        wb.name=paste(label, std, sep="."), wb.label=label))
  } else if (std=="medical"){
    return(new_waveband(w.low=290, w.high=320,
                        wb.name=paste(label, std, sep="."), wb.label=label))
  } else if (std=="none"){
    return(new_waveband(w.low=280, w.high=320,
                        wb.name=paste(label, std, sep="."), wb.label=label))
  } else {
    warning("Unsupported value '", std, "' supplied for 'std'.")
    return(NA)
  }
}
