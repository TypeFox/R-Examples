### R code of package 'pbivnorm', version 0.5-1 (2012-10-31).
### Author: Brenton Kenkel
###
### Based on code from 'mnormt' package by Adelchi Azzalini (University of
### Padua, Italy); used under GPL

##' Calculate probabilities from the CDF of a standard bivariate normal
##' distribution.
##'
##' This function returns values identical to those of \code{biv.nt.prob} in the
##' \pkg{mnormt} package, but is vectorized to reduce the number of Fortran
##' calls required for computation of many probabilities.
##' @title Standard bivariate normal CDF
##' @param x vector of upper integration limits for the CDF.  May also be a
##' two-column matrix, in which case \code{y} should not be used.
##' @param y vector of upper integration limits.
##' @param rho correlation parameter.
##' @param recycle whether to automatically recycle the vectors \code{x},
##' \code{y}, and \code{rho} to conform to whichever is longest.  If
##' \code{FALSE}, all three must be the same length.
##' @return Numeric vector of probabilities.
##' @references
##' Genz, A.  (1992).  Numerical Computation of Multivariate Normal
##' Probabilities.  \emph{J. Computational and Graphical Statist.}, \strong{1},
##' 141--149.
##'
##' Genz, A.  (1993).  Comparison of methods for the computation of multivariate
##' normal probabilities.  \emph{Computing Science and Statistics}, \strong{25},
##' 400--405.
##'
##' Genz, A.  Fortran code for \code{MVTDSTPACK} available at
##' \url{http://www.math.wsu.edu/math/faculty/genz/software/fort77/mvtdstpack.f}
##' (as of 2011-02-21).
##' @author Fortran code by Alan Genz (see references).  R interface by Brenton
##' Kenkel (\email{brenton.kenkel@@gmail.com}), based on code from Adelchi
##' Azzalini's \pkg{mnormt} package.
##' @examples
##' x <- rnorm(10)
##' y <- rnorm(10)
##' rho <- runif(10)
##'
##' pbivnorm(x, y, rho)
##'
##' X <- cbind(x, y)
##' pbivnorm(X, rho = rho)
##'
##' ## rho can be a single value, unless recycling is disallowed
##' rho <- runif(1)
##' pbivnorm(x, y, rho)
pbivnorm <- function(x, y, rho = 0, recycle = TRUE) {
    ## allow for x to be a two-column matrix
    if (length(dim(x))) {
        if (ncol(x) != 2)
            stop("'x' must have two columns if specified as a matrix")
        if (!missing(y) && !is.null(y))
            warning("'x' was specified as a matrix, so 'y' will be ignored")
        y <- x[, 2]
        x <- x[, 1]
    }

    ## sanity checks
    if (any(abs(rho) > 1))
        stop("'rho' must be a valid correlation (-1 <= rho <= 1)")

    ## turn infinite integration limits into their best finite equivalents
    x <- replace(x, x == Inf, .Machine$double.xmax)
    x <- replace(x, x == -Inf, -.Machine$double.xmax)
    y <- replace(y, y == Inf, .Machine$double.xmax)
    y <- replace(y, y == -Inf, -.Machine$double.xmax)

    ## Recycle so that `x`, `y`, and `rho` are the same length (if requested)
    lengths <- sapply(list(x, y, rho), length)
    if (recycle) {
        x <- rep(x, length.out = max(lengths))
        y <- rep(y, length.out = max(lengths))
        rho <- rep(rho, length.out = max(lengths))
    } else if (diff(range(lengths)) > 0) {
        stop("'x', 'y', and 'rho' must be the same length when recycling is disallowed")
    }

    ## coerce arguments to proper types to be passed to fortran
    lower <- as.double(c(0, 0))
    infin <- as.integer(c(0, 0))
    uppera <- as.double(x)
    upperb <- as.double(y)
    lt <- as.integer(length(x))
    prob <- double(lt)
    rho <- as.double(rho)

    ans <- .Fortran("PBIVNORM", prob, lower, uppera, upperb, infin, rho, lt,
                    PACKAGE="pbivnorm")[[1]]
    return(ans)
}
