################################################################################
### Part of the R package "polyCub".
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Copyright (C) 2009-2015 Sebastian Meyer
### Time-stamp: <[polyCub.exact.Gauss.R] 2015-02-15 17:27 (CET) by SM>
################################################################################


#' Quasi-Exact Cubature of the Bivariate Normal Density
#'
#' Integration is based on triangulation of the (transformed) polygonal domain
#' and formulae from the 
#' Abramowitz and Stegun (1972) handbook (Section 26.9, Example 9, pp. 956f.).
#' This method is quite cumbersome because the A&S formula is only for triangles
#' where one vertex is the origin (0,0). For each triangle of the
#' \code{\link[gpclib]{tristrip}} we have to check in which of the 6 outer 
#' regions of the triangle the origin (0,0) lies and adapt the signs in the 
#' formula appropriately. (AOB+BOC-AOC) or (AOB-AOC-BOC) or (AOB+AOC-BOC) or
#' (AOC+BOC-AOB) or \ldots. However, the most time consuming step is the
#' evaluation of \code{\link[mvtnorm]{pmvnorm}}.
#' 
#' @note The package \pkg{gpclib} (which is required to produce the
#' \code{tristrip}, since this is not yet implemented in \pkg{rgeos})
#' has a restricted license (commercial use prohibited).
#' It has to be accepted explicitly via
#' \code{\link{gpclibPermit}()} prior to using \code{polyCub.exact.Gauss}.
#'
#' @param polyregion a \code{"\link[rgeos:gpc.poly-class]{gpc.poly}"} polygon or
#' something that can be coerced to this class, e.g., an \code{"owin"} polygon
#' (converted via \code{\link{owin2gpc}} and -- given \pkg{rgeos} is available
#' -- \code{"SpatialPolygons"} also work. 
#' @param mean,Sigma mean and covariance matrix of the bivariate normal density
#' to be integrated.
#' @param plot logical indicating if an illustrative plot of the numerical
#' integration should be produced. Note that the \code{polyregion} will be
#' transformed (shifted and scaled).
#' @return The integral of the bivariate normal density over \code{polyregion}.
#' Two attributes are appended to the integral value:
#' \item{nEval}{
#' number of triangles over which the standard bivariate normal density had to 
#' be integrated, i.e. number of calls to \code{\link[mvtnorm]{pmvnorm}} and
#' \code{\link[stats]{pnorm}}, the former of which being the most time-consuming
#' operation.
#' }
#' \item{error}{
#' Approximate absolute integration error steming from the error introduced by
#' the \code{nEval} \code{\link[mvtnorm]{pmvnorm}} evaluations.
#' For this reason, the cubature method is in fact only
#' quasi-exact (as is the \code{pmvnorm} function).
#' }
#' @references
#' Abramowitz, M. and Stegun, I. A. (1972).
#' Handbook of Mathematical Functions with Formulas, Graphs, and Mathematical
#' Tables. New York: Dover Publications.
#' @keywords math spatial
#' @seealso \code{\link{circleCub.Gauss}} for quasi-exact cubature of the
#' isotropic Gaussian density over a circular domain.
#' @family polyCub-methods
#' @examples # see example(polyCub)
#' @import methods
#' @import sp
#' @importFrom stats cov2cor
#' @importFrom spatstat is.polygonal
#' @export
## NOTE: we don't import graphics::plot since it is already imported via sp

polyCub.exact.Gauss <- function (polyregion, mean = c(0,0), Sigma = diag(2),
                                 plot = FALSE)
{
    gpclibCheck(fatal=TRUE)
    if (is.polygonal(polyregion)) {
        polyregion <- owin2gpc(polyregion)
    } else if (!inherits(polyregion, "gpc.poly")) {
        if (inherits(polyregion, "SpatialPolygons") &&
            !requireNamespace("rgeos")) {
            stop("package ", sQuote("rgeos"), " is required to handle ",
                 "\"SpatialPolygons\" input")
        }
        polyregion <- as(polyregion, "gpc.poly")
    }
    
    ## coordinate transformation so that the standard bivariat normal density
    ## can be used in integrations (cf. formula 26.3.22)
    polyregion@pts <- transform_pts(polyregion@pts, mean = mean, Sigma = Sigma)
    
    ## triangulation: tristrip() returns a list where each element is a
    ## coordinate matrix of vertices of triangles 
    triangleSets <- gpclib::tristrip(polyregion)
    
### ILLUSTRATION ###
    if (plot) {
        plot(polyregion, poly.args=list(lwd=2), ann=FALSE)
        lapply(triangleSets, lines, lty=2)
    }
####################

    integrals <- vapply(X = triangleSets, FUN = function (triangles) {
        int <- 0
        error <- 0
        nTriangles <- nrow(triangles) - 2L
        for (i in seq_len(nTriangles)) {
            res <- .intTriangleAS(triangles[i+(0:2),])
            int <- int + res
            error <- error + attr(res, "error")
        }
        c(int, nTriangles, error)
    }, FUN.VALUE = numeric(3L), USE.NAMES = FALSE)
    int <- sum(integrals[1,])
    
    ## number of .V() evaluations (if there were no degenerate triangles)
    attr(int, "nEval") <- 6 * sum(integrals[2,])
    ## approximate absolute integration error
    attr(int, "error") <- sum(integrals[3,])
    return(int)
}



###########################
### Auxiliary Functions ###
###########################

## transform coordinates according to Formula 26.3.22
transform_pts <- function (pts, mean, Sigma)
{
    mx <- mean[1L]
    my <- mean[2L]
    rho <- cov2cor(Sigma)[1L,2L]
    sdx <- sqrt(Sigma[1L,1L])
    sdy <- sqrt(Sigma[2L,2L])
    lapply(pts, function (poly) {
        x0 <- (poly[["x"]] - mx) / sdx
        y0 <- (poly[["y"]] - my) / sdy
        list(x = (x0 + y0) / sqrt(2 + 2*rho),
             y = (y0 - x0) / sqrt(2 - 2*rho),
             hole = poly[["hole"]])
    })
}

## calculates the integral of the standard bivariat normal over a triangle ABC
.intTriangleAS <- function (xy)
{
    if (anyDuplicated(xy)) # degenerate triangle
        return(structure(0, error = 0))
    A <- xy[1,]
    B <- xy[2,]
    C <- xy[3,]
    intAOB <- .intTriangleAS0(A, B)
    intBOC <- .intTriangleAS0(B, C)
    intAOC <- .intTriangleAS0(A, C)
    
    # determine signs of integrals
    signAOB <- -1 + 2*.pointsOnSameSide(A,B,C)
    signBOC <- -1 + 2*.pointsOnSameSide(B,C,A)
    signAOC <- -1 + 2*.pointsOnSameSide(A,C,B)
    
    int <- signAOB*intAOB + signBOC*intBOC + signAOC*intAOC
    attr(int, "error") <- attr(intAOB, "error") +
        attr(intBOC, "error") + attr(intAOC, "error")
    return(int)
}

## calculates the integral of the standard bivariat normal over a triangle A0B
.intTriangleAS0 <- function (A, B)
{
    BmA <- B - A
    d <- sqrt(sum(BmA^2))
    h <- abs(B[2L]*A[1L] - A[2L]*B[1L]) / d   # distance of AB to the origin
    if (d == 0 || h == 0) # degenerate triangle: A == B or 0, A, B on a line
        return(structure(0, error = 0))
    
    k1 <- dotprod(A, BmA) / d
    k2 <- dotprod(B, BmA) / d
    V2 <- .V(h, abs(k2))
    V1 <- .V(h, abs(k1))
    res <- if (sign(k1) == sign(k2)) {
        ## A and B are on the same side of the normal line through 0
        abs(V2 - V1)
    } else {
        V2 + V1
    }
    attr(res, "error") <- attr(V1, "error") + attr(V2, "error")
    return(res)
}

## checks if point1 and point2 lie on the same side of a line through
## linepoint1 and linepoint2
.pointsOnSameSide <- function (linepoint1, linepoint2, point1, point2 = c(0,0))
{
    n <- c(-1,1) * rev(linepoint2-linepoint1)   # normal vector
    S <- dotprod(point1-linepoint1,n) * dotprod(point2-linepoint1,n)
    return(S > 0)
}

## calculates the integral of the standard bivariat normal
## over a triangle bounded by y=0, y=ax, x=h (cf. formula 26.3.23)
##' @importFrom stats pnorm
.V <- function(h,k) {
    if (k == 0) # degenerate triangle
        return(structure(0, error = 0))
    a <- k/h
    rho <- -a/sqrt(1+a^2)
    # V = 0.25 + L(h,0,rho) - L(0,0,rho) - Q(h) / 2
    # L(0,0,rho) = 0.25 + asin(rho) / (2*pi)
    # V = L(h,0,rho) - asin(rho)/(2*pi) - Q(h) / 2
    Lh0rho <- mvtnorm::pmvnorm(
        lower = c(h,0), upper = c(Inf,Inf),
        mean = c(0,0), corr = matrix(c(1,rho,rho,1), 2L, 2L)
    )
    Qh <- pnorm(h, mean = 0, sd = 1, lower.tail = FALSE)
    return(Lh0rho - asin(rho)/2/pi - Qh/2)
}
