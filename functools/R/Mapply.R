#' Wrapper for mapply function.
#'
#' Wrapper for \code{\link[base]{mapply}}.
#'
#' @param ... Arguments to vectorize over (vectors or lists of strictly positive length, or all of zero length).
#' @param .f A function to be applied.
#' @param more_args A list of other arguments to FUN.
#' @param simplify Logical or character string; attempt to reduce the result to a vector, matrix or higher dimensional array; see the simplify argument of \code{\link[base]{sapply}}.
#' @param use_names Logical; use names if the first ... argument has names, or if it is a character vector, use that character vector as the names.
#' @family functionals
#' @seealso \code{\link[base]{mapply}} for code and documentation.
#' @export
Mapply <- function(..., .f, more_args = NULL, simplify = TRUE, use_names = TRUE) {
  return(mapply(FUN = .f, ..., MoreArgs = more_args, SIMPLIFY = simplify, USE.NAMES = use_names))
}
