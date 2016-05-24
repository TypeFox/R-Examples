#' Definition of UV-A waveband
#'
#' UV-A according to CIE and ISO standrads: 315--400 nm.
#' UV-A according to common non-standard practice: 320--400 nm.
#'
#' @param std a character string "CIE", "ISO" or "none"
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' UVA()
#' UVA("none")
#' UVA("ISO")
#' UVA("CIE")
#'
#' @family unweighted wavebands
#'
UVA <- function(std="ISO") {
  label <- "UVA"
  if (std=="ISO" | std=="CIE"){
    return(new_waveband(w.low=315, w.high=400, wb.name=paste("UVA",std, sep="."), wb.label=label))
  } else if (std=="none"){
    return(new_waveband(w.low=320, w.high=400, wb.name=paste("UVA",std, sep="."), wb.label=label))
  } else {
    warning("Unsupported value for 'std'.")
    return(NA)
  }
}
