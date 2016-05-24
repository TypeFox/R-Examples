#' Working with \code{kdecopula} objects
#'
#' The function \code{\link{kdecop}} stores it's result in ojbect of class
#' \code{kdecopula}. The density estimate can be evaluated on arbitrary points
#' with \code{\link[kdecopula:dkdecop]{dkdecop}}; the cdf with
#' \code{\link[kdecopula:pkdecop]{pkdecop}}. Furthermore, synthetic data can be
#' simulated with \code{\link[kdecopula:rkdecop]{rkdecop}}.
#'
#' @aliases dkdecop pkdecop rkdecop
#'
#' @param u \code{mx2} matrix of evaluation points.
#' @param obj \code{kdecop} object.
#' @param stable logical; option for stabilizing the estimator: the estimated
#' density is cut off at \eqn{50}.
#' 
#' @return A numeric vector of the density/cdf or a \code{n x 2} matrix of 
#' simulated data.
#' 
#' @author Thomas Nagler
#' 
#' @seealso 
#' \code{\link[kdecopula:kdecop]{kdecop}},
#' \code{\link[kdecopula:plot.kdecopula]{plot.kdecopula}},
#' \code{\link[qrng:ghalton]{ghalton}}
#' 
#' @references 
#' Geenens, G., Charpentier, A., and Paindaveine, D. (2014).
#' Probit transformation for nonparametric kernel estimation of the copula
#' density.
#' arXiv:1404.4414 [stat.ME]. 
#' \cr \cr 
#' Nagler, T. (2014). 
#' Kernel Methods for Vine Copula Estimation.
#' Master's Thesis, Technische Universitaet Muenchen,
#' \url{https://mediatum.ub.tum.de/node?id=1231221} 
#' \cr \cr 
#' Cambou, T., Hofert, M., Lemieux, C. (2015).
#' A primer on quasi-random numbers for copula models, 
#' arXiv:1508.03483 [stat.CO]
#' 
#' 
#' @examples
#' 
#' ## load data and transform with empirical cdf
#' data(wdbc)
#' udat <- apply(wdbc[, -1], 2, function(x) rank(x)/(length(x)+1))
#' 
#' ## estimation of copula density of variables 5 and 6
#' dens.est <- kdecop(udat[, 5:6])
#' plot(dens.est) 
#' 
#' ## evaluate density estimate at (u1,u2)=(0.123,0.321)
#' dkdecop(c(0.123, 0.321), dens.est) 
#' 
#' ## evaluate cdf estimate at (u1,u2)=(0.123,0.321)
#' pkdecop(c(0.123, 0.321), dens.est) 
#' 
#' ## simulate 500 samples from density estimate
#' rkdecop(500, dens.est)
#'
dkdecop <- function(u, obj, stable = FALSE) {
    stopifnot(is.numeric(u))
    stopifnot(all(u > 0 & u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(is.logical(stable))

    ## define appropriately shaped udata frame for vapply
    u <- as.matrix(u)
    if (ncol(u) == 1)
        u <- matrix(u, 1L, nrow(u))

    ## adjust for flipping option of kdevine package
    d <- ncol(u)
    if (!is.null(obj$flip))
        u <- matrix(u[, 2:1], nrow(u), d)


    ## if independence copula is specified return 1
    if ("indep.copula" %in% class(obj))
        return(rep(1, nrow(u)))

    ## evaluate density  (use faster algorithm for d = 2)
    u <- pmin(pmax(u, 1e-3), 1 - 1e-3)
    if (d == 2) {
        out <- interp_2d(u,
                         obj$estimate,
                         obj$grid)
    } else {
        # define help indicators
        tmplst <- split(rep(seq(-1, 2, 1), d), ceiling(seq.int(4*d)/4))
        helpind <- as.matrix(do.call(expand.grid, tmplst))

        out <- interp(u,
                      obj$estimate,
                      obj$grid,
                      helpind)
    }

    ## stabilize output
    if (stable)
        out <- pmin(out, 10^(1 + d/2))

    ## return results
    out
}

#' @rdname dkdecop
#'
pkdecop <- function(u, obj) {
    stopifnot(all(u > 0 & u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    ## define appropriately shaped u matrix
    u <- as.matrix(u)
    if(ncol(u) == 1)
        u <- matrix(u, 1L, nrow(u))
    # adjust for flipping option of kdevine package
    if(!is.null(obj$flip)) {
        u <- matrix(u[, 2:1], nrow(u))
        cond.var <- ifelse(cond.var == 1, 2, 1)
    }
    d <- ncol(u)

    ## if independence copula is specified, return prod(u) directly
    if ("indep.copula" %in% class(obj))
        return(apply(u, 1, prod))

    ## define help objects
    tmplst <- split(rep(seq(-1, 2, 1), d), ceiling(seq.int(4*d)/4))
    helpind <- as.matrix(do.call(expand.grid, tmplst))
    m <- length(obj$grid)
    tmplst <- split(rep(obj$grid, d), ceiling(seq.int(m*d)/m))
    helpgrid <- as.matrix(do.call(expand.grid, tmplst))

    ## evaluate cdf
    eval_cdf(u,
             obj$estimate,
             obj$grid,
             helpgrid,
             helpind)
}

#' @param n integer; number of observations.
#' @param quasi logical; the default (\code{FALSE}) returns pseudo-random
#' numbers, use \code{TRUE} for quasi-random numbers (generalized Halton, see
#' \code{\link[qrng:ghalton]{ghalton}}).
#'
#' @rdname dkdecop
#'
rkdecop <- function(n, obj, quasi = FALSE) {
    n <- round(n)
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(is.logical(quasi))

    if (!quasi) {
        # simulate independent uniform random variables
        W <- cbind(runif(n), runif(n))
    } else {
        # generate quasi random numbers
        W <- ghalton(n, d = 2)
    }

    # if independence copula is specified, return W
    if("indep.copula" %in% class(obj))
        return(W)

    # invert h-function otherwise
    U2 <- inv_hfunc(W,
                    1L,
                    obj$estimate,
                    obj$grid)

    ## return results
    out <- cbind(W[, 1], U2)
    colnames(out) <- NULL
    out
}

# Evaluate h-function of kdecopula object (internal for kdevine-package)
#
# @param cond.var integer; either 1 or 2; specifies on which argument to
#  condition
#
hkdecop <- function(u, obj, cond.var) {
    stopifnot(all(u > 0) & all(u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(all(cond.var %in% c(1, 2)))
    ## define appropriately shaped u matrix
    u <- as.matrix(u)
    if(ncol(u) == 1)
        u <- matrix(u, 1L, nrow(u))
    # adjust for flipping option of kdevine package
    if(!is.null(obj$flip)) {
        u <- matrix(u[, 2:1], nrow(u))
        cond.var <- ifelse(cond.var == 1, 2, 1)
    }
    d <- ncol(u)

    ## if independence copula is specified, return the conditioned variable
    if ("indep.copula" %in% class(obj))
        return(u[, -cond.var])

    ## evaluate h-function otherwise
    if (d == 2) {
        # faster algorithm for d = 2
        return(eval_hfunc_2d(u,
                             as.integer(cond.var),
                             obj$estimate,
                             obj$grid))
    } else {
        # define help objects
        tmplst <- split(rep(seq(-1, 2, 1), d), ceiling(seq.int(4*d)/4))
        helpind <- as.matrix(do.call(expand.grid, tmplst))
        tmplst <- lapply(1:d,
                         prep_hfunc,
                         cond.var = cond.var,
                         grid = obj$grid,
                         d = d)
        helpgrid <- as.matrix(do.call(expand.grid, tmplst))
        uncond.var <- seq.int(d)[-cond.var]

        # call routine
        return(eval_hfunc(u,
                          as.integer(cond.var),
                          as.integer(uncond.var),
                          obj$estimate,
                          obj$grid,
                          helpgrid,
                          helpind))
    }
}

## Prepare the list for constructing the helpgrid in hkdecop
prep_hfunc <- function(i, cond.var, grid, d) {
    if (i %in% cond.var) {
        return(0)
    } else {
        return(grid)
    }
}

## Calculate inverse of h-function
hinvkdecop <- function(u, obj, cond.var) {
    stopifnot(all(u > 0 & u < 1))
    stopifnot(inherits(obj, "kdecopula"))
    stopifnot(all(cond.var %in% c(1, 2)))
    ## simulate independent uniform random variables
    W <- u

    # if independence copula is specified, return W
    if("indep.copula" %in% class(obj))
        return(W)

    ## invert h-function otherwise
    U2 <- inv_hfunc(W,
                    cond.var,
                    obj$estimate,
                    obj$grid)

    ## return results
    out <- cbind(W[, 1], U2)
    colnames(out) <- NULL
    out
}

