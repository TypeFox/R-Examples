#' Definition of UV-C waveband
#'
#' UV-C according to CIE and ISO standrads: 100--280 nm. UV-c according to
#' common non-standard practice: 200--280 nm. UV-C according to medical or
#' dermatological non-standard practice, e.g. Diffey (1991): 200--290 nm.
#'
#' @param std a character string "CIE", "ISO", "none", or "medical".
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @export
#'
#' @seealso \code{\link{new_waveband}}  \code{\link{waveband}}
#'
#' @examples
#' UVC()
#' UVC("CIE")
#' UVC("ISO")
#' UVC("none")
#' UVC("medical")
#'
#' @family unweighted wavebands
#'
UVC <- function(std="ISO") {
  label <- "UVC"
  if (std=="ISO" || std=="CIE"){
    return(new_waveband(w.low=100, w.high=280, wb.name=paste(label ,std, sep="."), wb.label=label))
  } else if (std=="medical"){
    return(new_waveband(w.low=220, w.high=290, wb.name=paste(label,std, sep="."), wb.label=label))
  } else if (std=="none"){
    return(new_waveband(w.low=200, w.high=280, wb.name=paste(label,std, sep="."), wb.label=label))
  } else {
    warning("Unsupported value '", std, "' supplied for 'std'.")
    return(NA)
  }
}
