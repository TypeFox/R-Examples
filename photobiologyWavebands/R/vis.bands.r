#' Definition of list of VIS wavebands
#'
#' Defined according to "ISO".
#'
#' @param std a character string "ISO" (ignored)
#' @return a list of wavebands
#' @export
#'
#' @examples
#' VIS_bands()
#' VIS_bands("ISO")
#'
#' @family lists of unweighted wavebands
#'
VIS_bands <- function(std="ISO"){
  return(list(Purple(std), Blue(std), Green(std), Yellow(std), Orange(std), Red(std)))
}
