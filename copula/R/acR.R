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


### McNeil, Neslehova (2009): Distribution function of the radial part #########

##' \bar{F}_R(x) = \sum_{k=0}^{d-1} (x^k * (-1)^k * psi^{(k)}(x)) / k!
##'
##' @title Distribution Function of the Radial Part for Archimedean Copulas
##' @param x evaluation points
##' @param family Archimedean family
##' @param theta parameter
##' @param d dimension
##' @param lower.tail logical indicating whether the survival function of F_R is
##'        returned
##' @param log.p logical indicating whether the logarithm is returned
##' @param ... additional arguments passed to absdPsi()
##' @return F_R at x
##' @author Marius Hofert
pacR <- function(x, family, theta, d, lower.tail = TRUE, log.p = FALSE, ...)
{
    ## basic checks
    stopifnot(family %in% .ac.longNames, d >= 1, (n <- length(x)) >= 1, x >= 0)
    family <- match.arg(family, choices=.ac.longNames)

    ## d == 1
    cop <- onacopulaL(family, list(theta, 1:d))
    psi.x <- cop@copula@psi(x, theta=theta)
    if(d == 1)
        return( if(lower.tail) if(log.p) log1p(-psi.x) else 1-psi.x
                else if(log.p) log(psi.x) else psi.x )

    ## d >= 2; compute \log\bar{F}_R(x)
    lpsi <- log(psi.x) # n-vector; k==0
    k <- 1:(d-1)
    lx <- log(x)
    lkf <- lfactorial(k) # length d-1
    mat1d1 <- vapply(k, function(k.) # (n, d-1) matrix; k = 1,..,d-1
                     cop@copula@absdPsi(x, theta=theta, degree=k., log=TRUE, ...) +
                     k.*lx - lkf[k.], rep(NA_real_, n))
    lsummands <- matrix(c(lpsi, mat1d1), nrow=n, ncol=d) # (n, d) matrix; also works for n=1 (!)

    ## compute log(\bar{F}_R(x))
    lFb <- lsum(t(lsummands)) # length n
    is0 <- x == 0
    if(any(is0)) lFb[is0] <- 0

    ## return
    if(log.p) if(lower.tail) log1mexp(-lFb) else lFb
    else if(lower.tail) -expm1(lFb) else exp(lFb)
}

##' Computing the quantile function of pacR()
##'
##' Compute x = F_R^{-1}(p) by solving F_R(x) = p w.r.t. p
##' Equivalently, solve log(1-F_R(x)) = log1p(-p) w.r.t. x
##' (numerically more stable; what we do here)
##' For log.p = TRUE, compute x = log(F_R^{-1}(exp(p)))
##' @title Computing the Quantile Function of pacR()
##' @param p probability in (0,1) where to evaluate qacR()
##' @param family Archimedean family (e.g., "AMH", "Clayton", "Frank", "Gumbel", "Joe")
##' @param theta parameter
##' @param d dimension
##' @param log.p logical; if TRUE, probabilities p are given as log(p).
##' @param interval interval for uniroot() [theoretically [0, Inf))
##' @param tol see ?uniroot()
##' @param maxiter see ?uniroot()
##' @param ... additional arguments passed to pacR()
##' @return the quantile function of F_R at x
##' @author Marius Hofert
qacR <- function(p, family, theta, d, log.p = FALSE, interval,
                 tol=.Machine$double.eps^0.25, maxiter=1000, ...)
{
    stopifnot(family %in% .ac.longNames, d >= 1, (n <- length(p)) >= 1, 0 < p,  p < 1)
    if(log.p) p <- exp(p) # not very efficient yet
    y <- log1p(-p) # to be able to work with highest precision with pacR()
    f <- function(x, y) pacR(x, family=family, theta=theta, d=d, lower.tail=FALSE, log.p=TRUE, ...) - y
    vapply(y, function(y.) uniroot(f, interval=interval, y=y.,
                                   tol=tol, maxiter=maxiter)$root, NA_real_)
}
