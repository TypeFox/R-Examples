################################################################################
### Part of the R package "polyCub".
### Free software under the terms of the GNU General Public License, version 2,
### a copy of which is available at http://www.r-project.org/Licenses/.
###
### Copyright (C) 2013-2015 Sebastian Meyer
### Time-stamp: <[polyCub.iso.R] 2015-02-16 12:18 (CET) by SM>
################################################################################


#' Cubature of Isotropic Functions over Polygonal Domains
#'
#' Conducts numerical integration of a two-dimensional isotropic function
#' \eqn{f(x,y) = f_r(||(x,y)-\boldsymbol{\mu}||)}{f(x,y) = f_r(||(x,y)-\mu||)},
#' with \eqn{\mu} being the center of isotropy, over a polygonal domain. 
#' It internally solves a line integral along the polygon boundary using
#' \code{\link{integrate}} where the integrand requires the antiderivative of
#' \eqn{r f_r(r)}), which ideally is analytically available and supplied to the
#' function as argument \code{intrfr}.
#' The two-dimensional integration problem thereby reduces to an efficient
#' adaptive quadrature in one dimension.
#' See Meyer and Held (2014, Supplement B, Section 2.4) for mathematical
#' details.
#'
#' @inheritParams plotpolyf
#' @param intrfr analytical antiderivative of \eqn{r f_r(r)} from 0 to \code{R}
#' (first argument, not necessarily named \code{"R"}, must be vectorized).
#' If missing, \code{intrfr} is approximated numerically using
#' \code{\link{integrate}} configured with \code{control}.
#' @param ... further arguments for \code{f} or \code{intrfr}.
#' @param center numeric vector of length 2, the center of isotropy.
#' @param control list of arguments passed to \code{\link{integrate}}, the
#' quadrature rule used for the line integral along the polygon boundary.
#' @param check.intrfr logical (or numeric vector) indicating if
#' (for which \code{r}'s) the supplied \code{intrfr} function should be
#' checked against a numeric approximation. This check requires \code{f}
#' to be specified. If \code{TRUE}, the set of test
#' \code{r}'s defaults to a \code{\link{seq}} of length 20 from 1 to
#' the maximum absolute x or y coordinate of any edge of the \code{polyregion}.
#' @param plot logical indicating if an image of the function should be plotted
#' together with the polygonal domain, i.e.,
#' \code{\link{plotpolyf}(polyregion, f, \dots)}.
#' @return The approximate integral of the isotropic function
#' \code{f} over \code{polyregion}.\cr
#' If the \code{intrfr} function is provided (which is assumed to be exact), an
#' upperbound for the absolute integration error is appended to the result as
#' attribute \code{"abs.error"}. It equals the sum of the absolute errors
#' reported by all \code{\link{integrate}} calls
#' (there is one for each edge of \code{polyregion}).
#' @author Sebastian Meyer
#' 
#' The basic mathematical formulation of this efficient integration for radially
#' symmetric functions was ascertained with great support by
#' Emil Hedevang (2013), Dept. of Mathematics, Aarhus University, Denmark.
#' @references
#' Hedevang, E. (2013). Personal communication at the Summer School on Topics in
#' Space-Time Modeling and Inference (May 2013, Aalborg, Denmark).
#'
#' Meyer, S. and Held, L. (2014).
#' Power-law models for infectious disease spread.
#' \emph{The Annals of Applied Statistics}, \bold{8} (3), 1612-1639.\cr
#' DOI-Link: \url{http://dx.doi.org/10.1214/14-AOAS743},
#' \href{http://arxiv.org/abs/1308.5115}{arXiv:1308.5115}
#' @keywords math spatial
#' @family polyCub-methods
#' @example inst/examples/polyCub.iso.R
#' @importFrom stats integrate
#' @export

polyCub.iso <- function (polyregion, f, intrfr, ..., center,
                         control = list(), check.intrfr = FALSE, plot = FALSE)
{
    polys <- xylist(polyregion) # transform to something like "owin$bdry"
                                # which means anticlockwise vertex order with
                                # first vertex not repeated
    getError <- !missing(intrfr) # can't estimate error of double approximation
    center <- as.vector(center, mode = "numeric")
    stopifnot(length(center) == 2L, is.finite(center))

    ## check 'intrfr' function
    rs <- if (isTRUE(check.intrfr)) {
        seq(1, max(abs(unlist(lapply(polys, "[", c("x","y"))))), length.out=20L)
    } else if (identical(check.intrfr, FALSE)) {
        numeric(0L)
    } else {
        check.intrfr
    }
    intrfr <- checkintrfr(intrfr, f, ..., center=center, control=control, rs=rs)

    ## plot polygon and function image
    if (plot) plotpolyf(polys, f, ...)
    
    ## do the cubature over all polygons of the 'polys' list
    .polyCub.iso(polys, intrfr, ..., center=center,
                 control=control, .witherror=getError)
}


##' Check the Integral of \eqn{r f_r(r)}
##'
##' This function is auxiliary to \code{\link{polyCub.iso}}.
##' The (analytical) integral of \eqn{r f_r(r)} from 0 to \eqn{R} is checked
##' against a numeric approximation using \code{\link{integrate}} for various
##' values of the upper bound \eqn{R}. A warning is issued if inconsistencies
##' are found.
##'
##' @inheritParams polyCub.iso
##' @param rs numeric vector of upper bounds for which to check the validity of
##' \code{intrfr}. If it has length 0, no checks are performed.
##' @param tolerance of \code{\link{all.equal.numeric}} when comparing
##' \code{intrfr} results with numerical integration. Defaults to the
##' relative tolerance used for \code{integrate}.
##' @return The \code{intrfr} function. If it was not supplied, its quadrature
##' version using \code{integrate} is returned.
##' @importFrom stats integrate
##' @export
checkintrfr <- function (intrfr, f, ..., center, control = list(),
                         rs = numeric(0L), tolerance = control$rel.tol)
{
    doCheck <- length(rs) > 0L
    if (!missing(f)) {
        f <- match.fun(f)
        rfr <- function (r, ...)
            r * f(cbind(center[1L]+r, center[2L], deparse.level=0L), ...)
        quadrfr1 <- function (R, ...) integrate(rfr, 0, R, ...)$value
        if (length(control))
            body(quadrfr1)[[2L]] <- as.call(c(as.list(body(quadrfr1)[[2L]]),
                                             control))
        quadrfr <- function (R, ...)
            vapply(X = R, FUN = quadrfr1, FUN.VALUE = 0, ..., USE.NAMES = FALSE)
        if (missing(intrfr)) {
            return(quadrfr)
        } else if (doCheck) {
            cat("Checking 'intrfr' against a numeric approximation ... ")
            stopifnot(is.vector(rs, mode="numeric"))
            if (is.null(tolerance))
                tolerance <- eval(formals(integrate)$rel.tol)
            ana <- intrfr(rs, ...)
            num <- quadrfr(rs, ...)
            if (!isTRUE(comp <- all.equal(num, ana, tolerance=tolerance))) {
                cat("\n->", comp, "\n")
                warning("'intrfr' might be incorrect: ", comp)
            } else cat("OK\n")
        }
    } else if (doCheck) {
        stop("numerical verification of 'intrfr' requires 'f'")
    }
    
    match.fun(intrfr)
}


##' \code{.polyCub.iso} is a \dQuote{bare-bone} version of \code{polyCub.iso}.
##' @rdname polyCub.iso
##' @param polys something like \code{owin$bdry}, but see \code{\link{xylist}}.
##' @param .witherror logical indicating if an upperbound for the absolute
##' integration error should be attached as an attribute to the result?
##' @export
.polyCub.iso <- function (polys, intrfr, ..., center,
                          control = list(), .witherror = FALSE)
{
    ints <- lapply(polys, polyCub1.iso,
                   intrfr, ..., center=center,
                   control=control, .witherror=.witherror)
    if (.witherror) {
        res <- sum(vapply(X=ints, FUN="[", FUN.VALUE=0, 1L, USE.NAMES=FALSE))
        attr(res, "abs.error") <-
            sum(vapply(X=ints, FUN="[", FUN.VALUE=0, 2L, USE.NAMES=FALSE))
        res
    } else {
        sum(unlist(ints, recursive=FALSE, use.names=FALSE))
    }
}

## cubature method for a single polygon
polyCub1.iso <- function (poly, intrfr, ..., center,
                          control = list(), .witherror = TRUE)
{
    xy <- cbind(poly[["x"]], poly[["y"]], deparse.level=0L)
    nedges <- nrow(xy)
    intedges <- erredges <- numeric(nedges)
    for (i in seq_len(nedges)) {
        v0 <- xy[i, ] - center
        v1 <- xy[if (i==nedges) 1L else i+1L, ] - center
        int <- lineInt(v0, v1, intrfr, ..., control=control)
        intedges[i] <- int$value
        erredges[i] <- int$abs.error
    }
    int <- sum(intedges)
    ## if (!is.null(poly$hole) && !isTRUE(all.equal(0, int))) {
    ##     if ((1 - 2 * as.numeric(poly$hole)) * sign(int) == -1)
    ##         warning("wrong sign if positive integral")
    ## }
    if (.witherror) {
        c(int, sum(erredges))
    } else {
        int
    }
}

## line integral for one edge
##' @importFrom stats integrate
lineInt <- function (v0, v1, intrfr, ..., control = list())
{
    d <- v1 - v0
    num <- v1[2L]*v0[1L] - v1[1L]*v0[2L]  # = d[2]*p[,1] - d[1]*p[,2]
                                          # for any point p on the edge
    if (num == 0) { # i.e., if 'center' is part of this polygon edge
        return(list(value = 0, abs.error = 0))
    }
    integrand <- function (t) {
        ## get the points on the edge corresponding to t
        p <- cbind(v0[1L] + t*d[1L], v0[2L] + t*d[2L], deparse.level=0L)
        norm2 <- .rowSums(p^2, length(t), 2L)
        ints <- intrfr(sqrt(norm2), ...)
        ##ints[is.infinite(ints)] <- 1e300
        num * ints / norm2
    }
    if (length(control)) {              # use slower do.call()-construct
        do.call("integrate", c(list(integrand, 0, 1), control))
    } else {
        integrate(integrand, 0, 1)
    }
}

## equally fast method _only_ for convex polygonal domains including the origin
## (formula obtained via polar coordinate representation)
lineInt2 <- function (v0, v1, intrfr, ..., control = list())
{
    d <- v1 - v0
    ld <- vecnorm(d)
    l0 <- vecnorm(v0)
    l1 <- vecnorm(v1)
    dp <- dotprod(v0,v1)
    theta <- acos((l0 - dp/l0) / ld)
    num <- sin(theta) * l0
    phispan <- acos(dp / l0 / l1)
    integrand <- function (phi, ...) {
        r <- num / sin(theta+phi)
        intrfr(r, ...)
    }
    do.call("integrate", c(list(integrand, 0, phispan, ...), control))
}
