#' Wrapper for vapply function.
#'
#' Wrapper for \code{\link[base]{vapply}}.
#'
#' @param .x A vector.
#' @param .f A function to be applied.
#' @param fun_value A (generalized) vector; a template for the return value from .f.
#' @param ... Optional arguments to .f.
#' @param use_names Logical; if TRUE and if X is character, use .x as names for the result unless it had names already.
#' @family functionals
#' @seealso \code{\link[base]{vapply}} for code and documentation.
#' @export
Vapply <- function(.x, .f, fun_value, ..., use_names = TRUE) {
  return(vapply(X = .x, FUN = .f, FUN.VALUE = fun_value, ..., USE.NAMES = use_names))
}
