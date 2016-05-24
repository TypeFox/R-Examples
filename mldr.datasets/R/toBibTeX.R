#' BibTeX entry associated to an mldr object
#' @description Gets the content of the \code{bibtex} member of the \code{mldr} object and returns it
#' @param object The mldr object whose BibTeX entry is needed
#' @param ... Additional parameters from the generic toBibtex function not used by toBibtex.mldr
#' @return A string with the BibTeX entry
#' @examples
#' \dontrun{
#' library(mldr.datasets)
#' cat(toBibtex(emotions))
#' }
#' @import utils
#' @export
toBibtex.mldr <- function(object, ...) {
  object$bibtex
}
