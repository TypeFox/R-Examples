#' Wrapper for tapply function.
#'
#' Wrapper for \code{\link[base]{tapply}}.
#'
#' @param .f A function to be applied.
#' @param .x A vector.
#' @param index List of one or more factors, each of same length as .x. The elements are coerced to factors by as.factor.
#' @param ... Optional arguments to .f.
#' @param simplify If FALSE, tapply always returns an array of mode "list". If TRUE (the default), then if FUN always returns a scalar, tapply returns an array with the mode of the scalar.
#' @family functionals
#' @seealso \code{\link[base]{tapply}} for code and documentation.
#' @export
Tapply <- function(.x, index, .f = NULL, ..., simplify = TRUE) {
  return(tapply(X = .x, INDEX = index, FUN = .f, ..., simplify = simplify))
}
