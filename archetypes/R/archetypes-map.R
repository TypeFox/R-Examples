
#' Archetypal maps
#'
#' Two-dimensional projection of the observations based on the alpha
#' coefficients into a space spanned by the (projected) archetypes.
#'
#' @param object An \code{\link{archetypes}} object
#' @param projection Projection function; see
#'   \code{\link{archmap_projections}}
#' @param projection_args Arguments passed to the projection function;
#'   see \code{\link{archmap_projections}}
#' @param rotate Rotation angle to rotate the projection
#' @param cex Character expansion of archetypes
#' @param col Color of observations
#' @param pch Point character of observations
#' @param xlab A label for the x-axis
#' @param ylab A label for the y-axis
#' @param axes Logical value to draw axes or not
#' @param asp The y/x aspect ratio
#' @param ... Arguments passed to the underlying plot function
#'
#' @return
#'   Invisible matrix with the projected archetypes
#'
#' @examples
#' \dontrun{
#'   data("skel", package = "archetypes")
#'   skel2 <- subset(skel, select = -Gender)
#'
#'   set.seed(1981)
#'   a <- archetypes(skel2, k = 5)
#'
#'   ## Simplex projection:
#'   archmap(a, col = skel$Gender)
#'
#'   ## Simplex projection with archetypes arranged according to their
#'   ## distances:
#'   archmap(a, col = skel$Gender,
#'           projection = tspsimplex_projection)
#'   archmap(a, col = skel$Gender,
#'           projection = tspsimplex_projection,
#'           projection_args = list(equidist = TRUE))
#'
#'   ## MDS projection:
#'   archmap(a, col = skel$Gender,
#'           projection = atypes_projection)
#' }
#'
#' @family archmap
#'
#' @export
archmap <- function(object, projection = simplex_projection,
                    projection_args = list(), rotate = 0,
                    cex = 1.5, col = 1, pch = 1, xlab = "", ylab = "",
                    axes = FALSE, asp = TRUE, ...) {
    
    .Deprecated("simplexplot", old = "archmap")

    stopifnot("archetypes" %in% class(object))
    stopifnot(is.function(projection))

    k <- object$k
    if( k < 3) {
      stop("Need at least 3 archetypes.\n")
    }

    ## Projection:
    cmds <- do.call(projection, c(list(parameters(object)), projection_args))


    ## Rotation:
    if ( rotate != 0 ){
      a <- pi*rotate/180
      A <- matrix(c(cos(a), -sin(a), sin(a),
                    cos(a)), ncol=2)
      cmds <- cmds %*% A
    }

    ## Active archetypes:
    hmds <- chull(cmds)
    active <- 1:k %in% hmds

    ## Plot region spanned by the projected archetypes:
    plot(cmds, type = "n", xlab = xlab, ylab = ylab, axes = axes, asp = asp, ...)
    points(coef(object) %*% cmds, col = col, pch = pch)

    rad <- ceiling(log10(k)) + 1.5
    polygon(cmds[hmds,])
    points(cmds[active,], pch=21, cex=rad*cex, bg="grey")
    text(cmds[active,], labels=(1:k)[active], cex=cex)
    if(any(!active)){
      points(cmds[!active,,drop=FALSE], pch=21, cex=rad*cex,
             bg="white", fg="grey")
      text(cmds[!active,,drop=FALSE], labels=(1:k)[!active],
           cex=cex, col="grey20")
    }

    invisible(cmds)
}



#' Archetypal map projections
#'
#' @param x Archetypes matrix
#' @param r Radius of the simplex projection
#'
#' @return
#'   Matrix with the projected archetypes
#'
#' @family archmap
#'
#' @aliases archmap_projections
#' @rdname archmap_projections
#' @export
simplex_projection <- function(x, r = 10) {
  phi <- seq(-pi, pi, length.out = nrow(x) + 1)
  phi <- phi[-1]

  cbind(x = r * cos(phi),
        y = r * sin(phi))
}



#' @param equidist Arrange archetypes equidistantly or in relation to
#'   their distance
#' @param ... Parameters for the \code{\link[TSP]{solve_TSP}} function
#' @rdname archmap_projections
#' @export
tspsimplex_projection <- function(x, r = 10, equidist = FALSE, ...) {
  stopifnot(require("TSP"))

  d <- dist(x)
  xo <- as.integer(solve_TSP(TSP(d), ...))

  if ( equidist ) {
    phi <- seq(-pi, pi, length.out = nrow(x) + 1)
    phi <- phi[-1][xo]
  } else {
    d <- as.matrix(d)
    phi <- mapply(function(i, j) d[i, j], xo, c(tail(xo, -1), xo[1]))
    phi <- cumsum((phi / sum(phi)) * 360)
    phi <- c(0, head(phi, -1))
    phi <- ((phi * 2 * pi) / 360) - pi
  }

  cbind(x = r * cos(phi),
        y = r * sin(phi))
}



#' @rdname archmap_projections
#' @export
atypes_projection <- function(x) {
  cmdscale(dist(x))
}
