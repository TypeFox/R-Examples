#' Wrapper for sapply function.
#'
#' Wrapper for \code{\link[base]{sapply}}.
#'
#' @param .x A vector.
#' @param .f A function to be applied.
#' @param simplify Logical or character string; should the result be simplified to a vector, matrix or higher dimensional array if possible?
#' @param use_names Logical; if TRUE and if .x is character, use .x as names for the result unless it had names already.
#' @param ... Optional arguments to f.
#' @family functionals
#' @seealso \code{\link[base]{sapply}} for code and documentation.
#' @export
Sapply <- function(.x, .f, ..., simplify = TRUE, use_names = TRUE) {
  return(sapply(X = .x, FUN = .f, ..., simplify = simplify, USE.NAMES = use_names))
}
