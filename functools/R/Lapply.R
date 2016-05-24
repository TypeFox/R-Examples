#' Wrapper for lapply function.
#'
#' Wrapper for \code{\link[base]{lapply}}.
#'
#' @param .x A vector.
#' @param .f A function to be applied.
#' @param ... Optional arguments to .f.
#' @family functionals
#' @seealso \code{\link[base]{lapply}} for code and documentation.
#' @export
Lapply <- function(.x, .f, ...) {
  return(lapply(X = .x, FUN = .f, ...))
}
