################################################################################
### Part of the R package "polyCub".
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Copyright (C) 2009-2014 Sebastian Meyer
### Time-stamp: <[polyCub.SV.R] 2014-09-27 14:10 (CEST) by SM>
################################################################################


##' Product Gauss Cubature over Polygonal Domains
##'
##' Product Gauss cubature over polygons as proposed by
##' Sommariva and Vianello (2007).
##'
##' @inheritParams plotpolyf
##' @param f a two-dimensional real function (or \code{NULL} to only compute
##' nodes and weights).
##' As its first argument it must take a coordinate matrix, i.e., a
##' numeric matrix with two columns, and it must return a numeric vector of
##' length the number of coordinates.
##' @param nGQ degree of the one-dimensional Gauss-Legendre quadrature rule
##' (default: 20) as implemented in function \code{\link[statmod]{gauss.quad}}
##' of package \pkg{statmod}. Nodes and weights up to \code{nGQ=60} are cached
##' in \pkg{polyCub}, for larger degrees \pkg{statmod} is required.
##' @param alpha base-line of the (rotated) polygon at \eqn{x = \alpha} (see
##' Sommariva and Vianello (2007) for an explication). If \code{NULL} (default),
##' the midpoint of the x-range of each polygon is chosen if no \code{rotation}
##' is performed, and otherwise the \eqn{x}-coordinate of the rotated point
##' \code{"P"} (see \code{rotation}). If \code{f} has its maximum value at the
##' origin \eqn{(0,0)}, e.g., the bivariate Gaussian density with zero mean,
##' \code{alpha = 0} is a reasonable choice.
##' @param rotation logical (default: \code{FALSE}) or a list of points
##' \code{"P"} and \code{"Q"} describing the preferred direction. If
##' \code{TRUE}, the polygon is rotated according to the vertices \code{"P"} and
##' \code{"Q"}, which are farthest apart (see Sommariva and Vianello, 2007). For
##' convex polygons, this rotation guarantees that all nodes fall inside the
##' polygon.
##' @param engine character string specifying the implementation to use.
##' Up to \pkg{polyCub} version 0.4-3, the two-dimensional nodes and weights 
##' were computed by \R functions and these are still available by setting
##' \code{engine = "R"}.
##' The new C-implementation is now the default (\code{engine = "C"}) and 
##' requires approximately 30\% less computation time.\cr
##' The special setting \code{engine = "C+reduce"} will discard redundant nodes
##' at (0,0) with zero weight resulting from edges on the base-line
##' \eqn{x = \alpha} or orthogonal to it. 
##' This extra cleaning is only worth its cost for computationally intensive
##' functions \code{f} over polygons which really have some edges on the
##' baseline or parallel to the x-axis.  Note that the old \R
##' implementation does not have such unset zero nodes and weights.
##' @param plot logical indicating if an illustrative plot of the numerical
##' integration should be produced.
##' @return The approximated value of the integral of \code{f} over
##' \code{polyregion}.\cr
##' In the case \code{f = NULL}, only the computed nodes and weights are
##' returned in a list of length the number of polygons of \code{polyregion},
##' where each component is a list with \code{nodes} (a numeric matrix with
##' two columns), \code{weights} (a numeric vector of length
##' \code{nrow(nodes)}), the rotation \code{angle}, and \code{alpha}.
##' @author Sebastian Meyer\cr
##' The product Gauss cubature is based on the
##' original \acronym{MATLAB} implementation \code{polygauss} by Sommariva and
##' Vianello (2007), which is available under the GNU GPL (>=2) license from
##' \url{http://www.math.unipd.it/~alvise/software.html}.
##' @references
##' Sommariva, A. and Vianello, M. (2007).
##' Product Gauss cubature over polygons based on Green's integration formula.
##' \emph{Bit Numerical Mathematics}, \bold{47} (2), 441-453.
##' @keywords math spatial
##' @family polyCub-methods
##' @importFrom graphics points
##' @examples # see example(polyCub)
##' @export

polyCub.SV <- function (polyregion, f, ...,
                        nGQ = 20, alpha = NULL, rotation = FALSE, engine = "C",
                        plot = FALSE)
{
    polys <- xylist(polyregion) # transform to something like "owin$bdry"
                                # which means anticlockwise vertex order with
                                # first vertex not repeated
    stopifnot(isScalar(nGQ), nGQ > 0,
              is.null(alpha) || (isScalar(alpha) && !is.na(alpha)))

    ## COMPUTE NODES AND WEIGHTS OF 1D GAUSS QUADRATURE RULE.
    ## DEGREE "N" (as requested) (ORDER GAUSS PRIMITIVE)
    nw_N <- gauss.quad(nGQ)
    ## DEGREE "M" = N+1 (ORDER GAUSS INTEGRATION)
    nw_M <- gauss.quad(nGQ + 1)

    ## Special case f=NULL: compute and return nodes and weights only
    if (is.null(f)) {
        return(lapply(X = polys, FUN = polygauss, nw_MN = c(nw_M, nw_N),
                      alpha = alpha, rotation = rotation, engine = engine))
    }
    
    ## Cubature over every single polygon of the "polys" list
    f <- match.fun(f)
    int1 <- function (poly) {
        nw <- polygauss(poly, c(nw_M, nw_N), alpha, rotation, engine)
        fvals <- f(nw$nodes, ...)
        cubature_val <- sum(nw$weights * fvals)
        ## if (!isTRUE(all.equal(0, cubature_val))) {
        ## if ((1 - 2 * as.numeric(poly$hole)) * sign(cubature_val) == -1)
        ## warning("wrong sign if positive integral")
        ## }
        cubature_val
    }
    respolys <- vapply(X=polys, FUN=int1, FUN.VALUE=0, USE.NAMES=FALSE)
    int <- sum(respolys)

### ILLUSTRATION ###
    if (plot) {
        plotpolyf(polys, f, ..., use.lattice=FALSE)
        for (i in seq_along(polys)) {
            nw <- polygauss(polys[[i]], c(nw_M, nw_N), alpha, rotation, engine)
            points(nw$nodes, cex=0.6, pch = i) #, col=1+(nw$weights<=0)
        }
    }
###################

    int
}

## this wrapper provides a partially memoized version of
## unname(statmod::gauss.quad(n, kind="legendre"))
gauss.quad <- function (n)
{
    if (n <= 61) { # results cached in R/sysdata.rda
        .NWGL[[n]]
    } else if (requireNamespace("statmod")) {
        unname(statmod::gauss.quad(n = n, kind = "legendre"))
    } else {
        stop("package ", sQuote("statmod"), " is required for nGQ > 60")
    }
}


##' Calculate 2D Nodes and Weights of the Product Gauss Cubature
##'
##' @param xy list with elements \code{"x"} and \code{"y"} containing the
##' polygon vertices in \emph{anticlockwise} order (otherwise the result of the
##' cubature will have a negative sign) with first vertex not repeated at the
##' end (like \code{owin.object$bdry}).
##' @param nw_MN unnamed list of nodes and weights of one-dimensional Gauss
##' quadrature rules of degrees \eqn{N} and \eqn{M=N+1} (as returned by
##' \code{\link[statmod]{gauss.quad}}): \code{list(s_M, w_M, s_N, w_N)}.
##' @inheritParams polyCub.SV
##' @references
##' Sommariva, A. and Vianello, M. (2007):
##' Product Gauss cubature over polygons based on Green's integration formula.
##' \emph{Bit Numerical Mathematics}, \bold{47} (2), 441-453.
##' @keywords internal
##' @useDynLib polyCub C_polygauss

polygauss <- function (xy, nw_MN, alpha = NULL, rotation = FALSE, engine = "C")
{
    ## POLYGON ROTATION
    
    xyrot <- if (identical(FALSE, rotation)) {
        if (is.null(alpha)) { # choose midpoint of x-range
            xrange <- range(xy[["x"]])
            alpha <- (xrange[1L] + xrange[2L]) / 2
        }
        angle <- 0
        xy[c("x", "y")]
    } else {
        ## convert to coordinate matrix
        xy <- cbind(xy[["x"]], xy[["y"]], deparse.level=0)
        ## determine P and Q
        if (identical(TRUE, rotation)) { # automatic choice of rotation angle
            ## such that for a convex polygon all nodes fall inside the polygon
            QP <- vertexpairmaxdist(xy)
            Q <- QP[1L,,drop=TRUE]
            P <- QP[2L,,drop=TRUE]
        } else if (is.list(rotation)) {  # predefined rotation
            P <- rotation$P
            Q <- rotation$Q
            stopifnot(is.vector(P, mode="numeric") && length(P) == 2L,
                      is.vector(Q, mode="numeric") && length(Q) == 2L)
            stopifnot(any(P != Q))
            rotation <- TRUE
        } else {
            stop("'rotation' must be logical or a list of points \"P\" and \"Q\"")
        }
        rotmat <- rotmatPQ(P,Q)
        angle <- attr(rotmat, "angle")
        if (is.null(alpha)) {
            Prot <- rotmat %*% P
            alpha <- Prot[1]
        }
        xyrot <- xy %*% t(rotmat)   # = t(rotmat %*% t(xy))
        ## convert back to list
        list(x = xyrot[,1L,drop=TRUE], y = xyrot[,2L,drop=TRUE])
    }

    ## number of vertices
    L <- length(xyrot[[1L]])
    
    
    ## COMPUTE 2D NODES AND WEIGHTS.
    
    if (engine == "R") {
        
        toIdx <- c(seq.int(2, L), 1L)
        nwlist <- mapply(.polygauss.side,
                         xyrot[[1L]], xyrot[[2L]],
                         xyrot[[1L]][toIdx], xyrot[[2L]][toIdx],
                         MoreArgs = c(nw_MN, alpha),
                         SIMPLIFY = FALSE, USE.NAMES = FALSE)

        nodes <- c(lapply(nwlist, "[[", 1L),
                   lapply(nwlist, "[[", 2L),
                   recursive=TRUE)
        dim(nodes) <- c(length(nodes)/2, 2L)
        weights <- unlist(lapply(nwlist, "[[", 3L),
                          recursive=FALSE, use.names=FALSE)

    } else { # use C-implementation

        ## degrees of cubature and vector template for results
        M <- length(nw_MN[[1L]])
        N <- length(nw_MN[[3L]])
        zerovec <- double(L*M*N)

        ## rock'n'roll
        nwlist <- .C(C_polygauss,
                     as.double(xyrot[[1L]]), as.double(xyrot[[2L]]),
                     as.double(nw_MN[[1L]]), as.double(nw_MN[[2L]]),
                     as.double(nw_MN[[3L]]), as.double(nw_MN[[4L]]),
                     as.double(alpha),
                     as.integer(L), as.integer(M), as.integer(N),
                     x = zerovec, y = zerovec, w = zerovec)[c("x", "y", "w")]
        
        nodes <- cbind(nwlist[[1L]], nwlist[[2L]], deparse.level=0)
        weights <- nwlist[[3L]]

        ## remove unset nodes from edges on baseline or orthogonal to it
        ## (note that the R implementation does not return such redundant nodes)
        if (engine == "C+reduce" && any(unset <- weights == 0)) {
            nodes <- nodes[!unset,]
            weights <- weights[!unset]
        }

    }

    ## back-transform rotated nodes by t(t(rotmat) %*% t(nodes))
    ## (inverse of rotation matrix is its transpose)
    list(nodes = if (rotation) nodes %*% rotmat else nodes,
         weights = weights, angle = angle, alpha = alpha)
}


## The working horse .polygauss.side below is an R translation
## of the original MATLAB implementation by Sommariva and Vianello (2007).

.polygauss.side <- function (x1, y1, x2, y2, s_loc, w_loc, s_N, w_N, alpha)
{
    if ((x1 == alpha && x2 == alpha) || (y2 == y1))
        ## side lies on base-line or is orthogonal to it -> skip
        return(NULL)
    
    if (x2 == x1) { # side is parallel to base-line => degree N
        s_loc <- s_N
        w_loc <- w_N
    }
    
    half_pt_x <- (x1+x2)/2
    half_length_x <- (x2-x1)/2
    
    half_pt_y <- (y1+y2)/2
    half_length_y <- (y2-y1)/2
    
    ## GAUSSIAN POINTS ON THE SIDE.
    x_gauss_side <- half_pt_x + half_length_x * s_loc
    y_gauss_side <- half_pt_y + half_length_y * s_loc

    scaling_fact_minus <- (x_gauss_side - alpha) / 2

    ## construct nodes and weights: x and y coordinates ARE STORED IN MATRICES.
    ## A COUPLE WITH THE SAME INDEX IS A POINT, i.e. P_i=(x(k),y(k)).
    ## Return in an unnamed list of nodes_x, nodes_y, weights
    ## (there is no need for c(nodes_x) and c(weights))
    list(alpha + tcrossprod(scaling_fact_minus, s_N + 1), # degree_loc x N
         rep.int(y_gauss_side, length(s_N)),              # length: degree_loc*N
         tcrossprod(half_length_y*scaling_fact_minus*w_loc, w_N)) # degree_loc x N
}

## NOTE: The above .polygauss.side() function is already efficient R code.
##       Passing via C only at this deep level (see below) turned out to be
##       slower than staying with R! However, stepping into C already for
##       looping over the edges in polygauss() improves the speed.
## ## @useDynLib polyCub C_polygauss_side
## .polygauss.side <- function (x1, y1, x2, y2, s_M, w_M, s_N, w_N, alpha)
## {
##     if ((x1 == alpha && x2 == alpha) || (y2 == y1))
##         ## side lies on base-line or is orthogonal to it -> skip
##         return(NULL)
##
##     parallel2baseline <- x2 == x1  # side is parallel to base-line => degree N
##     M <- length(s_M)
##     N <- length(s_N)
##     loc <- if (parallel2baseline) N else M
##     zerovec <- double(loc * N)
##     .C(C_polygauss_side,
##        as.double(x1), as.double(y1), as.double(x2), as.double(y2),
##        as.double(if (parallel2baseline) s_N else s_M),
##        as.double(if (parallel2baseline) w_N else w_M),
##        as.double(s_N), as.double(w_N), as.double(alpha),
##        as.integer(loc), as.integer(N),
##        x = zerovec, y = zerovec, w = zerovec)[c("x", "y", "w")]
## }


##' @importFrom stats dist
vertexpairmaxdist <- function (xy)
{
    ## compute euclidean distance matrix
    distances <- dist(xy)
    size <- attr(distances, "Size")
    
    ## select two points with maximum distance
    maxdistidx <- which.max(distances)
    lowertri <- seq_along(distances) == maxdistidx
    mat <- matrix(FALSE, size, size)
    mat[lower.tri(mat)] <- lowertri
    QPidx <- which(mat, arr.ind=TRUE, useNames=FALSE)[1L,]
    xy[QPidx,]    
}

rotmatPQ <- function (P, Q)
{
    direction_axis <- (Q-P) / sqrt(sum((Q-P)^2))
    
    ## determine rotation angle
    rot_angle_x <- acos(direction_axis[1L])
    rot_angle_y <- acos(direction_axis[2L])
    
    rot_angle <- if (rot_angle_y <= pi/2) {
        if (rot_angle_x <= pi/2) -rot_angle_y else rot_angle_y
    } else {
        if (rot_angle_x <= pi/2) pi-rot_angle_y else rot_angle_y
    }
    ## cat(sprintf(' [ANGLE CLOCKWISE (IN DEGREES)]: %5.5f\n', rot_angle*180/pi))

    ## rotation matrix
    rot_matrix <- diag(cos(rot_angle), nrow=2L)
    rot_matrix[2:3] <- c(-1,1) * sin(rot_angle) # clockwise rotation
    structure(rot_matrix, angle=rot_angle)
}
