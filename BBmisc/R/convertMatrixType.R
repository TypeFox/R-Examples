#' Converts storage type of a matrix.
#'
#' Works by setting \code{\link{mode}}.
#'
#' @param x [\code{matrix}]\cr.
#'   Matrix to convert.
#' @param type [\code{character(1)}]\cr
#'   New storage type.
#' @return [\code{matrix}].
#' @note \code{as.mytype} drops dimension when used on a matrix.
#' @export
convertMatrixType = function(x, type) {
  assertMatrix(x)
  assertChoice(type, c("integer", "numeric", "complex", "character", "logical"))
  storage.mode(x) = type
  return(x)
}
