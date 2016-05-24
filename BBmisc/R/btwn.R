#' Check if some values are covered by the range of the values in a second vector.
#'
#' @param x [\code{numeric(n)}]\cr
#'   Value(s) that should be within the range of \code{y}.
#' @param y [\code{numeric}]\cr
#'   Numeric vector which defines the range.
#' @return [\code{logical(n)}]. For each value in \code{x}: Is it in the range of \code{y}?
#' @usage x \%btwn\% y
#' @rdname btwn
#' @examples
#' x = 3
#' y = c(-1,2,5)
#' x %btwn% y
#' @export
`%btwn%` = function(x, y) {
  r = range(y)
  x <= r[2] & x >= r[1]
}
