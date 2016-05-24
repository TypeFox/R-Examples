#' @include pcplot.R
{}


#' Parallel coordinates of data and archetypes.
#' @param x An \code{\link{archetypes}} object.
#' @param data A matrix or data frame.
#' @param data.col Color of data lines.
#' @param data.lwd Width of data lines.
#' @param atypes.col Color of archetypes lines.
#' @param atypes.lwd Width of archetypes lines.
#' @param atypes.lty Type of archetypes lines.
#' @param chull An integer vector giving the indices of the points from
#'   \code{data} lying on the convex hull.
#' @param chull.col Color of convex hull lines.
#' @param chull.lwd Width of convex hull lines.
#' @param chull.lty Type of convex hull lines.
#' @param ... Passed to \code{\link{pcplot}} and \code{\link{lines.pcplot}}.
#' @return Undefined.
#' @method pcplot archetypes
#' @S3method pcplot archetypes
pcplot.archetypes <- function(x, data, data.col=gray(0.7), data.lwd=1,
                              atypes.col=2, atypes.lwd=2, atypes.lty=1,
                              chull=NULL, chull.col=1, chull.lwd=2, chull.lty=1, ...) {

  pcplot(data, col=data.col, lwd=data.lwd, ...)

  if ( !is.null(chull) )
    lines.pcplot(data[chull,], data,
                 col=chull.col, lwd=chull.lwd, lty=chull.lty, ...)

  lines.pcplot(parameters(x), data,
               col=atypes.col, lwd=atypes.lwd, lty=atypes.lty, ...)
}

