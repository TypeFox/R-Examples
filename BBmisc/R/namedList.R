#' Create named list, possibly initialized with a certain element.
#'
#' @param names [\code{character}]\cr
#'   Names of elements.
#' @param init [valid R expression]\cr
#'   If given all list elements are initialized to this, otherwise
#'   \code{NULL} is used.
#' @return [\code{list}].
#' @export
#' @examples
#' namedList(c("a", "b"))
#' namedList(c("a", "b"), init = 1)
namedList = function(names, init) {
  if (missing(names))
    return(list())
  n = length(names)
  if (missing(init))
    xs = vector("list", n)
  else
    xs = replicate(n, init, simplify = FALSE)
  setNames(xs, names)
}
