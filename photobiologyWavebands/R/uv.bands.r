#' Definition of list of UV wavebands
#'
#' Defined according to "ISO".
#'
#' @param std a character string "ISO" or "none"
#' @return a list of wavebands
#' @export
#'
#' @seealso \code{\link[photobiology]{waveband}}
#'
#' @examples
#' UV_bands()
#' UV_bands("ISO")
#'
#' @family lists of unweighted wavebands
#'
UV_bands <- function(std="ISO"){
  if (std=="none") {
    stdc <- "ISO"
  } else {
    stdc <- std
  }
  return(list(UVC(stdc), UVB(std), UVA(std)))
}
