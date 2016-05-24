#' Definition of UV waveband
#'
#' UV: 100--400 nm.
#'
#' @param std "ISO" or "CIE"
#'
#' @return a waveband object wavelength defining a wavelength range.
#'
#' @references
#' ISO and CIE standards
#'
#' @export
#'
#' @seealso \code{\link{new_waveband}}  \code{\link{waveband}}
#'
#' @examples
#' UV()
#' UV("ISO")
#'
#' @family unweighted wavebands
#'
UV <- function(std="ISO") {
  new_waveband(w.low=100, w.high=400, wb.name=paste("UV", std, sep="."), wb.label="UV")
}
