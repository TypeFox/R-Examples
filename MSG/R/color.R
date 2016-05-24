#' Generate colors from a vector
#'
#' This functions generates a color vector from an input vector, which can be of
#' the class numeric or factor.
#' @rdname vec2col
#' @param vec the numeric or factor vector
#' @param n the number of colors to be generated from the palette
#' @param name the name of the palette
#' @return a vector of colors corresponding to the input vector
#' @export
#' @author Yihui Xie <\url{http://yihui.name}>
#' @examples
#' ## convert factor to colors
#' with(iris, plot(Petal.Length, Petal.Width, col = vec2col(Species), pch = 19))
#'
#' # another palette
#' with(iris, plot(Petal.Length, Petal.Width, col = vec2col(Species, name = 'Dark2'),
#' pch = 19))
#'
#' ## turn numeric values to colors
#' with(iris, plot(Petal.Length, Petal.Width, col = vec2col(Petal.Width), pch = 19))
vec2col = function(vec, n, name) {
  UseMethod('vec2col')
}
#' @rdname vec2col
#' @export
#' @importFrom RColorBrewer brewer.pal
vec2col.default = function(vec, n, name) {
  if (missing(n)) n = 9
  if (missing(name)) name = 'Blues'
  idx = findInterval(vec, seq(min(vec, na.rm = TRUE), max(vec, na.rm = TRUE),
                              length.out = n))
  brewer.pal(n, name)[idx]
}
#' @rdname vec2col
#' @export
vec2col.factor = function(vec, n, name) {
  if (missing(n)) n = length(levels(vec))
  if (missing(name)) name = 'Accent'
  brewer.pal(n, name)[as.integer(vec)]
}
