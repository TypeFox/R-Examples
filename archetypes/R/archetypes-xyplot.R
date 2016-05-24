

#' Two-dimensional plot.
#' @param x An object.
#' @param ... Further arguments.
#' @return Undefined.
#' @export
xyplot <- function(x, ...) {
  UseMethod('xyplot')
}



#' Helper function to calculate the approximated convex hull.
#' @param zs An \code{archetypes} object.
#' @return Matrix with the points.
#' @noRd
ahull <- function(zs) {
  a <- rbind(parameters(zs), parameters(zs)[1,])
  xc <- a[,1]; xm <- mean(xc)
  yc <- a[,2]; ym <- mean(yc)

  real <- xc - xm
  imag <- yc - ym
  angle <- atan2(imag, real)

  index <- order(angle)

  return(a[c(index, index[1]),])
}



#' Plot of two-dimensional data and archetypes.
#' @param x An \code{\link{archetypes}} object.
#' @param y A matrix or data frame.
#' @param data.col Color of data points.
#' @param data.pch Type of data points.
#' @param data.bg Background of data points.
#' @param atypes.col Color of archetypes points.
#' @param atypes.pch Type of archetypes points.
#' @param ahull.show Show approximated convex hull.
#' @param ahull.col Color of approximated convex hull line.
#' @param chull An integer vector giving the indices of the points from
#'   \code{data} lying on the convex hull.
#' @param chull.col Color of convex hull points.
#' @param chull.pch Type of convex hull points.
#' @param adata.show Show approximated data with link to the original
#'   data.
#' @param adata.col Color of approximated data points.
#' @param adata.pch Type of approximated data points.
#' @param link.col Color of link between approximated and original data
#'   points.
#' @param link.lty Line type of link between approximated and original
#'    data points.
#' @param ... Passed to the underlying plot functions.
#' @return Undefined.
#' @note The link between approximated and original data is based on an
#'   idea and Matlab source code of Bernard Pailthorpe.
#' @method xyplot archetypes
#' @export
xyplot.archetypes <- function(x, y,
                              data.col = 1, data.pch = 19, data.bg = NULL,
                              atypes.col = 2, atypes.pch = 19,
                              ahull.show = TRUE, ahull.col = atypes.col,
                              chull = NULL, chull.col = gray(0.7), chull.pch = 19,
                              adata.show = FALSE, adata.col = 3, adata.pch = 13,
                              link.col = data.col, link.lty = 1, ...) {

  zs <- x; data <- y;

  plot(data, col = data.col, pch = data.pch, bg = data.bg, ...)
  points(parameters(zs), col = atypes.col, pch = atypes.pch, ...)

  if ( !is.null(chull) ) {
    points(data[chull,], col = chull.col, pch = chull.pch, ...)
    lines(data[c(chull, chull[1]),], col = chull.col, ...)
  }

  if ( ahull.show )
    lines(ahull(zs), col = ahull.col)


  if ( adata.show ) {
    ### Based on an idea of Bernard Pailthorpe.
    adata <- fitted(zs)
    link.col <- rep(link.col, length = nrow(adata))
    link.lty <- rep(link.lty, length = nrow(adata))

    points(adata, col = adata.col, pch = adata.pch, ...)
    for ( i in seq_len(nrow(data)) )
      lines(rbind(data[i,], adata[i,]), col = link.col[i],
            lty = link.lty[i], ...)
  }

  invisible(NULL)
}



#' Plot of two-dimensional data and weighted archetypes.
#' @param x An \code{\link{archetypes}} object.
#' @param y A matrix or data frame.
#' @param data.col Color of data points.
#' @param data.pch Type of data points.
#' @param data.bg Background of data points.
#' @param link.col Color of link between approximated and original data
#'   points.
#' @param link.lty Line type of link between approximated and original
#'    data points.
#' @param weights.type Weights to display; see \code{\link{weights.archetypes}}.
#' @param ... Arguments of \code{\link{xyplot.archetypes}}.
#' @method xyplot weightedArchetypes
#' @S3method xyplot weightedArchetypes
xyplot.weightedArchetypes <- function(x, y, data.col = 1,
                                      data.pch = 21, data.bg = gray,
                                      link.col = NULL, link.lty = NULL,
                                      weights.type = 'weights', ...) {

  w <- weights(x, type = weights.type)
  if ( is.matrix(w) )
    w <- diag(w)
  w <- 1 - w

  if ( is.null(link.col) )
    link.col <- ifelse(w == 1, 1, data.bg(w))

  if ( is.null(link.lty) )
    link.lty <- ifelse(w == 1, 2, 1)

  if ( is.function(data.col) )
    data.col <- data.col(w)

  xyplot.archetypes(x, y, data.pch = data.pch,
                    data.col = data.col, data.bg = data.bg(w),
                    link.col = link.col, link.lty = link.lty, ...)
}



#' Plot of two-dimensional data and robust archetypes.
#' @param x An \code{\link{archetypes}} object.
#' @param y A matrix or data frame.
#' @param ... Arguments of \code{\link{xyplot.weightedArchetypes}} and
#'   \code{\link{xyplot.robustArchetypes}}
#' @method xyplot robustArchetypes
#' @S3method xyplot robustArchetypes
xyplot.robustArchetypes <- function(x, y, ...) {
  xyplot.weightedArchetypes(x, y, weights.type = 'reweights', ...)
}



#' Plot of two-dimensional data and stepArchetypes.
#' @param x An \code{\link{stepArchetypes}} object.
#' @param y A matrix or data frame.
#' @param data.col Color of data points.
#' @param data.pch Type of data points.
#' @param atypes.col Color of archetypes points.
#' @param atypes.pch Type of archetypes points.
#' @param ahull.show Show approximated convex hull.
#' @param ahull.col Color of approximated convex hull line.
#' @param ... Passed to the underlying plot functions.
#' @return Undefined.
#' @method xyplot stepArchetypes
#' @S3method xyplot stepArchetypes
#' @export
xyplot.stepArchetypes <- function(x, y,
                                  data.col=gray(0.7), data.pch=19,
                                  atypes.col=(seq_len(length(x) * length(x[[1]]))+1),
                                  atypes.pch=19, ahull.show=TRUE, ahull.col=atypes.col, ...) {

  zs <- x; data <- y;

  flatzs <- unlist(zs, recursive=FALSE)

  plot(data, col=data.col, pch=data.pch, ...)
  for ( i in seq_along(flatzs) ) {
    a <- flatzs[[i]]
    points(parameters(a), col=atypes.col[i], pch=atypes.pch, ...)

    if ( ahull.show )
      lines(ahull(a), col=ahull.col[i])
  }

  invisible(NULL)
}


