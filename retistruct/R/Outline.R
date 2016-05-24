##' Construct an outline object. This sanitises the input points
##' \code{P}, as described below.
##'
##' @title Outline constructor
##' @param P The points of the outline. The last point is not repeated.
##' @param scale The length of one unit of \code{P} in
##' micrometres. When images are present, this is the length of the
##' side of a pixel in the image.
##' @param im An image as a \code{raster} object
##' @return An \code{Outline} object containing the following:
##' \item{\code{P}}{A N-by-2 matrix of points of the \code{Outline} arranged in anticlockwise order}
##' \item{\code{gf}}{For each row of \code{P}, the index of \code{P} that is next in the outline travelling anticlockwise (forwards)}
##' \item{\code{gb}}{For each row of \code{P}, the index of \code{P} that is next in the outline travelling clockwise (backwards)}
##' \item{\code{im}}{The image as a \code{raster} object}
##' \item{\code{scale}}{The length of one unit of \code{P} in micrometres}
##' @author David Sterratt
Outline <- function(P, scale=NA, im=NULL) {
  o <- list()
  o$P <- P
  o$h <- 1:nrow(P)
  t <- TriangulatedOutline(o, n=NA)
  class(o) <- "outline"
  o$P <- t$P
  o$gf <- t$gf
  o$gb <- t$gb
  o$im <- im
  o$scale <- scale
  return(o)
}

##' Plot flat \code{\link{Outline}}. 
##'
##' @title Flat plot of outline
##' @param x \code{\link{Outline}} object
##' @param axt whether to plot axes
##' @param ylim y-limits
##' @param image If \code{TRUE} the image (if it is present) is
##' displayed behind the outline
##' @param scalebar If  numeric and if the Outline has a \code{scale}
##' field, a scale bar of length \code{scalebar} mm is plotted.  If
##' \code{scalebar} is \code{FALSE} or there is no scale information
##' in the \code{\link{Outline}} \code{x}  the scale bar is supressed.
##' @param add If \code{TRUE}, don't draw axes; add to existing plot.
##' @param lwd.outline Line width of outline
##' @param ... Other plotting parameters
##' @method flatplot outline
##' @author David Sterratt
##' @export
flatplot.outline <- function(x, axt="n", ylim=NULL,
                             image=TRUE,
                             scalebar=1,
                             add=FALSE,
                             lwd.outline=1,
                             ...) {
  plot.image <- image
  ## If there is no scale information, don't present a scale bar
  scalebar <- ifelse(is.numeric(scalebar) && !is.null(x$scale), scalebar, FALSE)

  with(x, {
    s <- which(!is.na(gb))                # source index
    d <- na.omit(gb)                      # destination index

    if (!add) {
      if (plot.image && !is.null(x$im)) {
        xs <- 1:ncol(im)
        ys <- 1:nrow(im)
        plot(NA, NA, xlim=c(0, max(xs)), ylim=c(0, max(ys)), asp=1,
             xaxt=axt, yaxt=axt, bty="n",
             xlab="", ylab="")
        ## rasterImage crashes on some systems, but not others.
        rasterImage(im, 0, 0, ncol(im), nrow(im))
      } else {
        xs <- P[s,1]
        ys <- P[s,2]
        plot(xs, ys, asp=1,
             pch=".", xaxt=axt, yaxt=axt, xlab="", ylab="",
             bty="n", ylim=ylim)
      }
    }
    segments(P[s,1], P[s,2], P[d,1], P[d,2],
             col=getOption("outline.col"), lwd=lwd.outline)

    ## Plot scalebar if required. scalebar is length in mm.
    if (!add && scalebar && !is.na(scale)) {
      sby <- min(ys) - 0.02*(max(ys) - min(ys))
      sblen <- 1000*scalebar/(scale)
      lines(c(max(xs) - sblen, max(xs)),c(sby, sby), lwd=2)
    }
  })
}

##' Simplify an outline object by removing verticies bordering short
##' edges while not encroaching on any of the outline. At present,
##' this is done by finding concave vertices. It is safe to remove
##' these, at the expense of increasing the area a bit.
##'
##' @title Simplify an outline object by removing short edges
##' @param o \code{outline} object to simplify
##' @param min.frac.length the minumum length as a fraction of the
##' total length of the outline. 
##' @param plot whether to display plotting or not during simplification
##' @return Simlified \code{outline} object
##' @author David Sterratt
simplify.outline <- function(o, min.frac.length=0.001, plot=FALSE) {
  P <- o$P
  N <- nrow(P)                        # Number of vertices
  Q <- rbind(P, P[1,])                # Convenience variable
  v <- diff(Q)                         # Vector of each edge
  l <- vecnorm(v)                     # Length of each edge
  ## Compute outer products at each vertex
  e <- extprod3d(cbind(v[c(N, 1:(N-1)),], 0), cbind(v, 0))[,3]
  
  ## Find short edges
  S <- l/sum(l) < min.frac.length

  ## Find indicies of points that can be removed.
  ## They have to be concave or colinear (e<=0). And they need to border short edges
  i.rem <- which((e <= 0) & (S | (S[c(N, 1:(N-1))])))

  if (plot) {
    ## Plot short edges...
    plot(P, col="white")
    if (any(S)) {
      segments(P[S,1], P[S,2], P[S,1]+v[S,1], P[S,2] + v[S,2], col="red")
    }
    ## and longer ones.
    if (any(!S)) {
      segments(P[!S,1], P[!S,2], P[!S,1]+v[!S,1], P[!S,2] + v[!S,2], col="black")
    }
    ## Plot colinear, convex and concave points
    points(P[e>0,1], P[e>0, 2], col="green")
    points(P[e==0,1], P[e==0, 2], col="orange")
    points(P[e<0,1], P[e<0, 2], col="blue")
    ## Plot points to remove
    points(P[i.rem,1], P[i.rem, 2], pch="X", col="brown")
  }

  ## If there are any points to remove, remove the first one
  if (length(i.rem) > 0) {
    message(paste("simplify.outline: Removing vertex", i.rem[1]))
    return(simplify.outline(list(P=P[-i.rem[1],], scale=o$scale, im=o$im)))
  } else {
    return(Outline(P, o$scale, o$im))
  }
}
