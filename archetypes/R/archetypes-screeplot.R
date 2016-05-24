

#' Screeplot of stepArchetypes.
#'
#' Screeplot draws the residual sum of square curve based on the best
#' model of each step.
#'
#' @param x A \code{\link{stepArchetypes}} object.
#' @param type Draw lines or a barplot.
#' @param ... Passed to underlying plot functions.
#' @return Undefined.
#' @importFrom stats screeplot
#' @method screeplot stepArchetypes
#' @S3method screeplot stepArchetypes
screeplot.stepArchetypes <- function(x, type=c('lines', 'barplot'), ...) {
  zs <- bestModel(x)

  a <- sapply(zs, nparameters)
  b <- sapply(zs, rss)

  if ( type[1] == 'lines' ) {
    plot(a, b, type='b', xlab='Archetypes', ylab='RSS', ...)
    axis(1, at=a, ...)
  }
  else {
    barplot(b, names.arg=a, xlab='Archetypes', ylab='RSS', ...)
  }
}
