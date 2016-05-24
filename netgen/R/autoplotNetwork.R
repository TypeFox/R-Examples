#' Autoplot function.
#'
#' Generates a \code{\link[ggplot2]{ggplot}} object. Nice possibility to
#' visualize 2-dimensional (clustered) networks in the euclidean plane.
#'
#' @param object [\code{Network}]\cr
#'   Network.
#' @param path [\code{integer}]\cr
#'   An integer vector containing the order of cities of a path. Keep in mind,
#'   that instances with \eqn{n} nodes and \eqn{m} depots have \eqn{n + m}
#'   coordinates, with the \eqn{1,\ldots,m} first coordinates belonging to
#'   the depots.
#' @param close.path [\code{logical(1)}]\cr
#'   Logical indicating whether the path passed by \code{path} should be
#'   closed to a cycle. Default is \code{FALSE}.
#' @param path.colour [\code{character(1)}]\cr
#'   Colour of the lines linking nodes on a path. Default is \dQuote{gray}.
#' @param ... [any]\cr
#'   Currently not used.
#' @return [\code{\link[ggplot2]{ggplot}}]
#' @examples
#' \dontrun{
#' # here we have no depots ...
#' x = generateClusteredNetwork(n.points = 30L, n.cluster = 2L)
#' pl = autoplot(x, path = 1:3)
#' # ... and here we have two depots: the path visits the depots in this case
#' x = generateRandomNetwork(n.points = 30L, n.depots = 2L)
#' pl = autoplot(x, path = 1:3, path.colour = "tomato")
#' }
#' @export
autoplot.Network = function(object,
  path = NULL, close.path = FALSE, path.colour = "gray",
  ...) {
  if (ncol(object$coordinates) > 2L) {
    stopf("Only 2-dimensional networks can be plotted.")
  }

  df = as.data.frame(object, include.extras = TRUE)

  if (testClass(object, "ClusteredNetwork")) {
    df$membership = as.factor(df$membership)
  }

  if (hasDepots(object)) {
    depot.idx = which(df$types == "depot")
    df.depots = df[depot.idx, , drop = FALSE]
  }

  pl = ggplot(data = df, mapping = aes_string(x = "x1", y = "x2"))

  if (!is.null(path)) {
    assertInteger(path, min.len = 2L, any.missing = FALSE)
    assertCharacter(path.colour, len = 1L, any.missing = FALSE)
    assertFlag(close.path)
    if (close.path) {
      path = c(path, path[1])
    }
    path.coords = df[path, , drop = FALSE]
    pl = pl + geom_path(data = path.coords, colour = path.colour)
  }

  if (!is.null(df$membership)) {
    pl = pl + geom_point(aes_string(colour = "membership"))
  } else {
    pl = pl + geom_point(colour = "tomato")
  }

  if (hasDepots(object)) {
    pl = pl + geom_point(data = df.depots, colour = "black", size = 4)
    pl = pl + geom_point(data = df.depots, colour = "white", size = 3)
  }
  pl = pl + ggtitle(as.character(object))
  pl = decorateGGPlot(pl, lower = object$lower, upper = object$upper)
  return(pl)
}
