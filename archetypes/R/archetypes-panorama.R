


#' Panorma plot for archetypes.
#' @param object An \code{\link{archetypes}}-related object.
#' @param data A matrix or data frame.
#' @param distfn Distance function.
#' @param xlab Label of xaxis.
#' @param ylab Label of yaxis.
#' @param order Order the distances.
#' @param col Color of distances.
#' @param pch Plot character of distances.
#' @param cex magnification of the distances.
#' @param atypes.col Color of archetype distances.
#' @param atypes.pch Plot character of archetype distances.
#' @param atypes.cex Magnification of the archetype distances.
#' @param ylim The y limits of the plot.
#' @param ... Passed to the underlying \code{plot} call.
#' @S3method panorama archetypes
#' @method panorama archetypes
#' @examples
#'   \dontrun{
#'   data(toy)
#'   a <- archetypes(toy, 3)
#'   panorama(a, toy)
#'
#'   ## See demo(robust-ozone).
#'   }
panorama.archetypes <- function(object, data, distfn = distEuclidean,
                                xlab = 'Index', ylab = 'Distance',
                                order = TRUE, col = 1, pch = 1, cex = 1,
                                atypes.col = (seq(length = nparameters(object)) + 1),
                                atypes.pch = rep(19, nparameters(object)),
                                atypes.cex = rep(1, nparameters(object)),
                                ylim = NULL, ...) {

  n1 <- nrow(data)
  n2 <- nparameters(object)

  data <- rbind(data, parameters(object))
  dist <- distfn(data, parameters(object))

  x <- seq(length = n1 + n2)

  col <- c(rep(col, n1), atypes.col)
  pch <- c(rep(pch, n1), atypes.pch)
  cex <- c(rep(cex, n1), atypes.cex)

  #if ( is.null(ref.order) )
  #  ix <- x
  #else
  #  ix <- order(dist[, ref.order])

  ix <- x
  r <- x

  opar <- par(no.readonly = TRUE)
  on.exit(par(opar))

  if ( is.null(ylim) )
    ylim <- c(0, max(dist))

  mar <- opar$mar
  mar[2] <- max(mar[2], 5)

  layout(matrix(seq(length = n2), nrow = n2, byrow = TRUE))
  par(mar = mar, cex = opar$cex)

  for ( i in seq(length = n2) ) {
    if ( order ) {
      ix <- order(dist[, i])
      r <- rank(dist[, i])
    }

    plot(x, dist[ix, i], ylim = ylim, xlab = xlab,
         ylab = ylab, col = col[ix], pch = pch[ix], cex = cex[ix], ...)

    or <- tail(r, n2)
    points(or, dist[tail(x, n2), i], pch = atypes.pch,
           cex = atypes.cex, col = atypes.col)

    mtext(sprintf('Archetype %s', i), side = 2, line = 4,
          cex = par('cex'))
  }


  invisible(dist)
}



#' Euclidean distance function (copied from flexclust)
#' @param x Data matrix.
#' @param centers Archetypes
#' @return Matrix with euclidean distance between each
#'   data point and each center.
#' @noRd
distEuclidean <- function (x, centers) {
    if (ncol(x) != ncol(centers))
        stop(sQuote("x"), " and ", sQuote("centers"), " must have the same number of columns")
    z <- matrix(0, nrow = nrow(x), ncol = nrow(centers))
    for (k in 1:nrow(centers)) {
        z[, k] <- sqrt(colSums((t(x) - centers[k, ])^2))
    }
    z
}
