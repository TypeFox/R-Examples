#' Wrapper for apply function.
#'
#' Wrapper for \code{\link[base]{apply}}.
#'
#' @param .x An array, including a matrix.
#' @param .f A function to be applied.
#' @param .margin A vector giving the subscripts which the function will be applied over.
#' @param ... Optional arguments to f.
#' @family functionals
#' @seealso \code{\link[base]{apply}} for code and documentation.
#' @export
Apply <- function(.x, .f, .margin, ...) {
  return(apply(X = .x, MARGIN = .margin, FUN = .f, ...))
}
