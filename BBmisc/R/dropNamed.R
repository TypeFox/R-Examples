#' Drop named elements of an object.
#'
#' @param x [any]\cr
#'   Object to drop named elements from.
#'   For a matrix or a data frames this function drops named columns via
#'   the second argument of the binary index operator \code{[,]}.
#'   Otherwise, the unary index operator \code{[]} is used for dropping.
#' @param drop [\code{character}]\cr
#'   Names of elements to drop.
#' @return Subset of object of same type as \code{x}. The object is not simplified,
#'   i.e, no dimensions are dropped as \code{[,,drop = FALSE]} is used.
#' @export
dropNamed = function(x, drop = character(0L)) {
  assertCharacter(drop, any.missing = FALSE)
  if (length(drop) == 0L)
    return(x)

  if (is.matrix(x) || is.data.frame(x))
    x[, setdiff(colnames(x), drop), drop = FALSE]
  else
    x[setdiff(names(x), drop)]
}
