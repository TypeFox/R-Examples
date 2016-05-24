## Copyright (C) 2012 Marius Hofert, Ivan Kojadinovic, Martin Maechler, and Jun Yan
##
## This program is free software; you can redistribute it and/or modify it under
## the terms of the GNU General Public License as published by the Free Software
## Foundation; either version 3 of the License, or (at your option) any later
## version.
##
## This program is distributed in the hope that it will be useful, but WITHOUT
## ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
## FOR A PARTICULAR PURPOSE. See the GNU General Public License for more
## details.
##
## You should have received a copy of the GNU General Public License along with
## this program; if not, see <http://www.gnu.org/licenses/>.


source(system.file("Rsource", "utils.R", package = "copula"))## ./utils.R
##--> tryCatch.W.E(), canGet()
stopifnot(requireNamespace("Runuran"))
## "import" what I need:
ud    <- Runuran::ud
ur    <- Runuran::ur
udgig <- Runuran::udgig
pinvd.new <- Runuran::pinvd.new

## stopifnot(require(optimx)) -- if needed load the Namespace only

### GIG generator and related functions ########################################

## generator (parameters theta_1 in IR, theta_2 in (0,Inf))
## note: theta_1 != 0 and theta_2 = 0 reduces to copClayton@psi(t, 1/theta_1)
psi.GIG <- function(t, theta)
    (1+t)^(-theta[1]/2) * besselK(theta[2]*sqrt(1+t),nu=theta[1])/
    besselK(theta[2],nu=theta[1])

## generator inverse
## note: the initial interval is only valid for theta_1 > -1/2
##       it can be large, but it's finite; it is based on an inequality about
##       modified Bessel functions of the third kind (besselK) given in
##       Paris (1984) ("An inequality for the Bessel function J_v(vx)",
##       SIAM J. Math. Anal. 15, 203--205)
iPsi.GIG <- function(t, theta, upper=NULL, ...) {
    if(is.null(upper)) {
	if(theta[1] > -0.5) upper <- function(x) (1-log(x)/theta[2])^2-1 else
        stop("initial interval for iPsi.GIG fails")
    }
    res <- numeric(length(t))
    is0 <- t == 0
    is1 <- t == 1
    n01 <- !(is0 | is1)
    res[is0] <- Inf
    res[is1] <- 0
    t. <- t[n01]
    up <- upper(t.)
    res[n01] <- unlist(lapply(seq_along(t.), function(i) {
	uniroot(function(t..) psi.GIG(t.., theta)-t.[i],
                interval=c(0, up[i]), ...)$root
    }))
    if(is.matrix(t)) matrix(res, ncol=ncol(t)) else res
}

## generator derivatives
absdPsi.GIG <- function(t, theta, degree=1, n.MC=0, log=FALSE) {
    res <- numeric(length(t))
    iInf <- is.infinite(t)
    res[iInf] <- -Inf # log(0)
    if(any(!iInf)) {
	t. <- t[!iInf]
        if(n.MC > 0) {
            V <- V0.GIG(n.MC, theta)
            res[!iInf] <- copula:::lsum(-V %*% t(t.) + degree*log(V) - log(n.MC))
        } else {
            res[!iInf] <- degree*log(theta[2]/2)-((theta[1]+degree)/2)*log1p(t.) +
                log(besselK(theta[2]*sqrt(1+t.), nu=theta[1]+degree, expon.scaled=TRUE)) -
                    log(besselK(theta[2], nu=theta[1], expon.scaled=TRUE)) -
                        (sqrt(1+t.)-1)*theta[2]
        }
    }
    r <- if(log) res else exp(res)
    if(is.matrix(t)) matrix(r, ncol=ncol(t)) else r
}

## absolute value of the derivative of the generator inverse -- NOWHERE USED !!
absdiPsi.GIG <- function(t, theta, log=FALSE) {
    absdPsi.GIG(iPsi.GIG(t, theta=theta), theta=theta, log=log)
}

## density of the GIG copula
dacopula.GIG <- function(u, theta, n.MC=0, log=FALSE) {
    if(!is.matrix(u)) u <- rbind(u)
    if((d <- ncol(u)) < 2) stop("u should be at least bivariate") # check that d >= 2
    ## f() := NaN outside and on the boundary of the unit hypercube
    res <- rep.int(NaN, n <- nrow(u))
    n01 <- apply(u,1,function(x) all(0 < x, x < 1)) # indices for which density has to be evaluated
    if(!any(n01)) return(res)
    ## auxiliary results
    u. <- u[n01,, drop=FALSE]
    psiI <- iPsi.GIG(u., theta=theta) # this is costly if d is large
    psiI.sum <- rowSums(psiI)
    ## main part
    if(n.MC > 0) { # Monte Carlo
        res[n01] <- absdPsi.GIG(psiI.sum, theta=theta, degree=d, n.MC=n.MC, log=TRUE) -
            rowSums(absdPsi.GIG(psiI, theta=theta, log=TRUE))
    } else { # explicit
        ## auxiliary function for density evaluation
        h <- function(t, k, theta, log=FALSE) {
            s1pt <- sqrt(1+t)
            B <- besselK(theta[2]*s1pt, nu=theta[1]+k)
            if(log) log(B)-(theta[1]+k)*log(s1pt) else B/s1pt^(theta[1]+k)
        }
        ## result
        res[n01] <- h(psiI.sum, k=d, theta=theta, log=TRUE) + (d-1)*
            log(besselK(theta[2], nu=theta[1])) - rowSums(h(psiI, k=1, theta=theta,
                                  log=TRUE))
    }
    if(log) res else exp(res)
}

## V0 random number generator
V0.GIG <- function(n, theta) {
    dens <- udgig(theta[1], 1, theta[2]^2) # the three args lambda, psi, chi for the GIG density as on page 497 in McNeil, Frey, Embrechts (2005)
    gen <- pinvd.new(dens) # works fast, via approximation of the quantile function by piecewise polynomials
    ur(gen, n)/2
}

## density of V0
dV0.GIG <- function(x, theta, log=FALSE) {
    dens <- udgig(theta[1], 1, theta[2]^2)
    if(log) log(2*ud(dens, 2*x)) else 2*ud(dens, 2*x)
}

## Kendall's tau for fixed theta_1 as a function in theta_2
## quiet == TRUE means that neither warnings nor errors are returned, just NAs
tau.GIG. <- function(theta, upper=c(200,200), quiet=TRUE, ...) {
    if(theta[1] < 0) stop("tau.GIG.: only implemented for nu >= 0")
    stopifnot(theta[1] <= upper[1], 0 < theta[2], theta[2] <= upper[2])
    ## determine minimal theta such that integration still works (before
    ## asymptotic is used)
    theta.min <- if(theta[1] < 50) {
	1e-4
    } else {
	if(theta[1] < 100) {
            0.1
        } else {
            l10 <- log(10)
            exp(l10/50*theta[1]-3*l10)
        }
    }
    ## either use heuristic or integration
    if(theta[2] < theta.min) {
        1/(1+2*theta[1])
    } else {
        tau.integrand <- function(t, theta) {
            st <- sqrt(1+t)
            t*(theta[2]*besselK(theta[2]*st, nu=theta[1]+1)/
               (st^(theta[1]+1)* besselK(theta[2], nu=theta[1])))^2
        }
        int <- tryCatch.W.E(integrate(tau.integrand, lower=0, upper=Inf, theta=theta,
                                      subdivisions=1000, ...)$value)
        warn.err <- inherits(int$value, what="simpleError") ||
            inherits(int$value, what="simpleWarning")
        ## are cautious: warning already gives NA:
        if(warn.err && quiet) NA else 1-int$value
    }
}

## wrapper for tau.GIG. [vectorized]
tau.GIG <- function(theta, upper=c(200,200), quiet=TRUE, ...) {
    if(!is.matrix(theta)) theta <- matrix(theta, ncol=2)
    apply(theta, 1L, tau.GIG., upper=upper, quiet=quiet, ...)
}

## inverse of Kendall's tau for all parameters fixed except the one given as NA
## note: initial interval has to be specified [e.g., c(1e-30,200)] since there
##       are no adequate bounds known
iTau.GIG <- function(tau, theta, upper=c(200,200), quiet=TRUE, iargs=list(), ...) {
    stopifnot(length(i <- which(is.na(theta))) == 1)
    tau.min <- tau.GIG(upper, upper=upper, quiet=quiet)
    if(tau <= tau.min) stop("iTau.GIG: tau must be greater than tau.min=", tau.min)
    interval <-
        if(i == 1) { # wanted: nu; given: theta
            c(0, upper[1])
        } else { # wanted: theta; given: nu
            if(theta[1] < 0) stop("iTau.GIG: only implemented for nu >= 0")
            tau.max <- 1/(1+2*theta[1])
            ## make sure we are in the range of attainable Kendall's tau
            if(tau >= tau.max)
                ## assumes tau to be falling in theta (numerical experiments show this behavior)
                stop("iTau.GIG: the supremum of attainable taus is ", round(tau.max, 8))
            c(1e-100, upper[2])
        }
    replace(theta, i, uniroot(function(th) {
        do.call(tau.GIG, c(list(replace(theta, i, th),
                                upper = upper,
                                quiet = quiet), iargs)) - tau },
        interval=interval, ...)$root)
}

## generate vectors of random variates from a GIG copula
rnacopula.GIG <- function(n, d, theta)
    psi.GIG(matrix(rexp(n*d), nrow=n, ncol=d)/V0.GIG(n, theta), theta)

