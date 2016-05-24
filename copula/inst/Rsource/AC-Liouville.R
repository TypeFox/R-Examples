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


## Note:
## - Random number generation (RNG) for simplex Archimedean copulas, Liouville
##   copulas, and Archimedean-Liouville copulas
## - Kendall's tau + inverse for simplex Archimedean copulas
## - Partly based on code from Alexander J. McNeil


### RNG for several ingredient distributions ###################################

##' @title Generating vectors of random variates from a d-dimensional simplex
##'        distribution
##' @param n sample size
##' @param d dimension
##' @return (n, d) matrix
##' @author Marius Hofert
rSimplex <- function(n, d) {
    stopifnot(n==as.integer(n), n>=0,
              d==as.integer(d), d>=0)
    E <- matrix(rexp(n*d), nrow=n, ncol=d)
    E/matrix(rowSums(E), nrow=n, ncol=d)
}

##' @title Generating vectors of random variates from a d-dimensional simplex
##'        distribution
##' @param n sample size
##' @param alpha vector of alphas for the Dirichlet distribution (positive integers)
##' @return (n, d) matrix
##' @author Marius Hofert
rDirichlet <- function(n, alpha) {
    stopifnot(n==as.integer(n), n>=0,
              alpha > 0)
    Z <- vapply(alpha, function(a) rgamma(n, a), numeric(n))
    Z/rowSums(Z)
}

##' Generating random variates from GPD(xi, beta) distribution
##'
##' The df of the GPD(xi, beta) distribution is
##'      G_{xi,beta}(x) = 1-(1+xi*x/beta)^(-1/xi) if xi!=0
##'                     = 1-exp(-x/beta)          if xi =0
##' for x>=0 when xi>=0 and x in [0,-beta/xi] when xi<0
##' @title Generating random variates from GPD(xi, beta) distribution
##' @param n sample size
##' @return n random variates
##' @author Marius Hofert
rGPD <- function(n, xi, beta) {
    stopifnot(n==as.integer(n), n>=1, length(xi)==1, length(beta)==1, beta>0)
    u <- runif(n)
    if(any(is0)) -beta*log1p(-u) else (beta/xi)*((1-u)^(-xi)-1)
}

##' @title Generating random variates from a Pareto distribution on [1,Inf)
##' @param n sample size
##' @return n random variates
##' @author Marius Hofert
rPareto <- function(n, kappa) (1-runif(n))^(-1/kappa)


### Kendall's tau + inverse for simplex Archimedean copulas ####################

##' @title Multivariate Kendall's tau for Simplex ACs
##' @param theta parameter
##' @param d dimension
##' @param Rdist distribution of the radial part
##' @param ... additional arguments passed to integrate()
##' @return Kendall's tau
##' @author Marius Hofert
##' @note *Multivariate* Kendall's tau according to Joe (1990)
tauACsimplex <- function(theta, d=2, Rdist=c("Gamma", "IGamma"), ...) {
    Rdist <- match.arg(Rdist)
    switch(Rdist,
           "Gamma" = {
               integrand <- function(y, d, th) (1-y)^(d-1) * y^(th-1) * (1+y)^(-2*th)
               ## note: numerically critical for small theta (large tau) or large theta
               Int <- integrate(integrand, lower=0, upper=1, d=d, th=theta)$value / beta(theta, theta)
               ## tau
               (2^d * Int - 1) / (2^(d-1) - 1)
           },
           "IGamma" =, "Pareto" =, "IPareto" = {
               stop("Not implemented yet; see McNeil, Neslehova (2010) for formulas")
           },
           stop("Not implemented yet"))
}

##' @title Inverse Multivariate Kendall's tau for Simplex Archimedean Copulas
##' @param tau Kendall's tau
##' @param d dimension
##' @param Rdist distribution of the radial part
##' @param interval theta-interval for uniroot()
##' @param ... additional arguments passed to uniroot()
##' @return theta such that tau(theta) = tau
##' @author Marius Hofert
##' @note non-critical numerical interval: c(1e-3, 1e2)
iTauACsimplex <- function(tau, d=2, Rdist=c("Gamma", "IGamma"), interval, ...)
    uniroot(function(th) tauACsimplex(th, d=d, Rdist=Rdist) - tau,
            interval=interval, ...)$root


### RNG for simplex Archimedean copulas ########################################

## auxiliary function for g():
## compute log(Gamma(a, x)), a in IR, x >= 0
## note: pgamma(x, a, lower=FALSE, log.p=TRUE) + lgamma(a) only works for a > 0
lg <- function(a, x) # vectorized in x
{
    stopifnot(x >= 0, length(a) == 1)
    n <- length(x)
    res <- numeric(n)
    i0 <- x==0
    if(any(i0)) res[i0]   <- lgamma(a)
    if(any(!i0)) res[!i0] <- log(gsl:::gamma_inc(a, x[!i0])) # if x>0 log(int_x^Inf t^(a-1)exp(-t)dt); no 'proper' log() exists
    res
}

## function from Example 9 of McNeil, Neslehova (2010) with k ~> d, n ~> k:
## g_{alpha, d, theta}(x) = sum_{k=1}^d (-x)^(d-k) \binom{d-1}{k-1} Gamma(k+theta-alpha, x)/Gamma(theta)
g <- function(x, alpha, d, theta){
    n <- length(x)
    k <- seq_len(d)
    sgns <- (-1)^(d-k) # d-vector
    ## remaining part: in log-scale
    lx <- outer(log(x), d-k) # (n, d)-matrix
    lc <- lchoose(d-1, k-1) # d-vector
    lgam <- sapply(k+theta-alpha, function(a) lg(a, x=x)) # (n, d)-matrix
    x.coeff <- exp(rep(lc, each=n) + lx + lgam - lgamma(theta)) # (n, d)-matrix
    rowSums(rep(sgns, each=n)*x.coeff)
}

##' @title Williamson d-transforms
##' @param t argument t>=0, a vector
##' @param d "dimension" d
##' @param theta parameter theta
##' @param Rdist distribution of the radial part
##' @return psi(t)
##' @author Marius Hofert
psiW <- function(t, d, theta, Rdist=c("Gamma", "IGamma", "Pareto", "IPareto"))
{
    stopifnot(t>=0, is.vector(t))
    Rdist <- match.arg(Rdist)
    switch(Rdist,
           "Gamma"={ # see Example 1 in McNeil, Neslehova (2010)
               n <- length(t)
               res <- numeric(n)
               i0 <- t==0
               if(any(i0))  res[i0]  <- 1
               if(any(!i0)) res[!i0] <- g(t[!i0], alpha=d, d=d, theta=theta)
               res
           },
           "IGamma"={ # see Example 2 in McNeil, Neslehova (2010)
               k <- 1:d
               arg <- d-k+theta # note: different from "Gamma"
               lg <- sapply(arg, function(a) pgamma(1/t, shape=a, log.p=TRUE) + lgamma(a)) # (length(t), length(k))-matrix
               lc <- lchoose(d-1, k-1)
               coeff <- exp(rep(lc, each=length(t)) + lg - lgamma(theta))
               x <- outer(-t, d-k, FUN="^") # (length(t), length(k))-matrix
               rowSums(coeff*x)
           },
           "Pareto"={ # see Example 3 in McNeil, Neslehova (2010)
               lth <- log(theta)
               lres <- lth - theta*log(t) + pbeta(pmin(1,x), shape1=theta, shape2=d, log.p=TRUE) +
                   lgamma(theta) + lgamma(d) - lgamma(theta+d)
               exp(lres)
           },
           "IPareto"={ # see Example 4 in McNeil, Neslehova (2010)
               S <- function(t, k, theta)
                   sapply(k, function(k.){
                       (-1)^(k.+1) * if(k.==theta) t^k.*log(t) else (t^theta-t^k.)/(theta-k.)
                   })
               k <- 1:d
               c. <- choose(d-1, k-1)
               S. <- S(t, k=k-1, theta=theta)
               if(!is.matrix(S.)) S. <- rbind(S., deparse.level=0L)
               theta * rowSums(rep(c., each=length(t)) * S.)
           },
           stop("wrong Rdist"))
}

##' @title Generating vectors of random variates from a d-dimensional simplex Archimedean
##'        distribution
##' @param n sample size
##' @param d dimension
##' @param theta parameter
##' @param Rdist distribution of the radial part
##' @return (n, d) matrix
##' @author Marius Hofert
rACsimplex <- function(n, d, theta, Rdist=c("Gamma", "IGamma", "Pareto", "IPareto"))
{
    stopifnot(n==as.integer(n), length(n)==1, n>=0,
              d==as.integer(d), length(d)==1, d>=0,
              length(theta)==1, theta>0) # might have to be adjusted if other families are added
    Rdist <- match.arg(Rdist)
    R <- switch(Rdist,
                "Gamma" = rgamma(n, theta),
                "IGamma" = 1/rgamma(n, theta),
                "Pareto" = rPareto(n, theta),
                "IPareto" = 1/rPareto(n, theta),
                stop("wrong Rdist"))
    S <- rSimplex(n, d=d)
    apply(R*S, 2, function(t) psiW(t, d=d, theta=theta, Rdist=Rdist))
}


### Kendall's tau + inverse for Liouville copulas ##############################

##' @title Kendall's tau for Liouville copulas
##' @param theta parameter (for Rdist)
##' @param alpha vector of alphas for the Dirichlet distribution (positive integers)
##' @param Rdist distribution of the radial part
##' @param n.MC Monte Carlo sample size
##' @return Kendall's tau
##' @author Marius Hofert
##' @note See Proposition 5 of McNeil, Neslehova (2010)
tauLiouville <- function(theta, alpha, Rdist="Gamma", n.MC=1e4) {
    stopifnot(alpha == as.integer(alpha), alpha >= 1, length(alpha) == 2)
    Rdist <- match.arg(Rdist)
    Y <- switch(Rdist,
                "Gamma" = {
                    R  <- rgamma(n.MC, theta)
                    R. <- rgamma(n.MC, theta)
                    R/R.
                },
                stop("Not implemented yet"))
    asum <- sum(alpha)
    ## EY1Y
    lEY1Y <- numeric(asum-1) # at least length 1
    for(l in seq_len(asum-1)) lEY1Y[l] <- log(mean(Y^(l-1) * pmax(1-Y, 0)^(asum-l)))
    ## compute the double sum (use index shift here)
    a1 <- alpha[1]
    a2 <- alpha[2]
    lga <- lgamma(asum)
    lb <- lbeta(a1, a2)
    S <- 0
    for(i in seq_len(a1)) {
        for(j in seq_len(a2)) {
            k <- i+j
            lfctr <- lbeta(a1+i-1, a2+j-1) + lga - lb - lfactorial(i-1) -
                lfactorial(j-1) - lgamma(asum-k+2)
            S <- S + exp(lfctr + lEY1Y[k-1]) # could be done more intelligently
        }
    }
    ## tau
    4*S-1
}

##' @title Inverse Kendall's tau for Liouville copulas
##' @param tau Kendall's tau
##' @param alpha vector of alphas for the Dirichlet distribution (positive integers)
##' @param Rdist distribution of the radial part
##' @param n.MC Monte Carlo sample size
##' @param interval theta-interval for uniroot()
##' @param tol tolerance for uniroot()
##' @param ... additional arguments passed to uniroot()
##' @return theta such that tau(theta) = tau
##' @author Marius Hofert
iTauLiouville <- function(tau, alpha, Rdist="Gamma", n.MC=1e4, interval, tol=1e-4, ...)
    uniroot(function(th) tauLiouville(th, alpha=alpha, Rdist=Rdist, n.MC=n.MC) - tau,
            interval=interval, tol=tol, ...)$root


### RNG for Liouville copulas (with given frailty distribution) ################

##' @title Marginal survival function of Liouville distributions
##' @param x numeric vector
##' @param j index of the marginal
##' @param alpha vector of alphas for the Dirichlet distribution (positive integers)
##' @param theta parameter theta
##' @param Rdist distribution of the radial part
##' @return vector of length as x
##' @author Marius Hofert
HbarL <- function(x, j, alpha, theta, Rdist=c("Gamma", "IGamma", "Clayton"))
{
    stopifnot(x>=0, alpha==as.integer(alpha), (d <- length(alpha))>=1,
              j==as.integer(j), length(j)==1, 1<=j, j<=d)
    Rdist <- match.arg(Rdist)
    switch(Rdist, # general formula: after Theorem 2 in McNeil, Neslehova (2010)
           "Gamma"={ # see Example 9 in McNeil, Neslehova (2010)
               ## bar(H)_j(x) = sum_{k=1}^{alpha_j} \binom{alpha-1}{k-1} x^{k-1} g_{alpha, alpha-k+1, theta}(x)
               asum <- sum(alpha)
               aj <- alpha[j]
               k <- seq_len(aj)
               c. <- choose(asum-1, k-1) # aj-vector
               x. <- outer(x, k-1, FUN="^") # (length(x), aj)-matrix
               g. <- sapply(asum-k+1, function(d.) g(x, alpha=asum, d=d., theta=theta)) # (length(x), aj)-matrix
               rowSums(rep(c., each=length(x))*x.*g.)
           },
           "IGamma"={ # see Example 10 in McNeil, Neslehova (2010)
               asum <- sum(alpha)
               aj <- alpha[j]
               k <- 1:aj
               c. <- choose(asum-1, k-1)
               x. <- outer(x, k-1, FUN="^") # (length(x), length(k))-matrix
               g. <- exp(lgamma(theta+k-1)-lgamma(theta))
               p. <- sapply(k-1, function(k.) psiW(x, d=asum-k., theta=theta+k., Rdist=Rdist))
               rowSums(rep(c.*g., each=length(x))*p.*x.)
           },
           "Clayton"={ # note: we use the simplified generator (!)
               aj <- alpha[j]
               k <- 0:(aj-1)
               lc. <- lgamma(1/theta+k)-lgamma(1/theta)
               coeff <- exp(rep(lc.-lfactorial(k), each=length(x))+outer(log(x/(1+x)), k))
               pmax(copClayton@psi(x, theta=theta), 0) * rowSums(coeff)
           },
           stop("wrong Rdist"))
}

##' @title Generating vectors of random variates from a Liouville copula
##' @param n sample size
##' @param alpha vector of alphas for the Dirichlet distribution (positive integers)
##' @param theta parameter theta
##' @param Rdist distribution of the radial part
##' @return (n, d=length(alpha)) matrix
##' @author Marius Hofert
##' @note See, for example, Example 9 and 10 in McNeil, Neslehova (2010)
rLiouville <- function(n, alpha, theta, Rdist=c("Gamma", "IGamma"))
{
    stopifnot(n==as.integer(n), length(n)==1, n>=0,
              (d <- length(alpha))>=1, length(theta)==1,
              theta>0) # might have to be adjusted if other families are added
    Rdist <- match.arg(Rdist)
    R <- switch(Rdist,
                "Gamma"=rgamma(n, theta),
                "IGamma"=1/rgamma(n, theta),
                stop("wrong Rdist"))
    D <- rDirichlet(n, alpha=alpha) # (n,d)
    X <- R*D # (n,d)
    sapply(1:d, function(j) HbarL(X[,j], j=j, alpha=alpha, theta=theta,
                                  Rdist=Rdist))
}


### RNG for Archimedean-Liouville copulas ######################################
### (with frailty distribution given by Williamson trafos) #####################

##' @title Generating vectors of random variates from an Archimedean-Liouville copula
##' @param n sample size
##' @param alpha vector of alphas for the Dirichlet distribution (positive integers)
##' @param theta parameter theta
##' @param family Archimedean family
##' @return (n, d=length(alpha)) matrix
##' @author Marius Hofert
##' @note See, for example, Example 11 in McNeil, Neslehova (2010)
rACLiouville <- function(n, alpha, theta, family=c("Clayton")){
    stopifnot(theta>=-1/((asum <- sum(alpha))-1),
              alpha==as.integer(alpha))
    if(theta <= 0) stop("Negative theta not yet implemented")
    d <- length(alpha)
    family <- match.arg(family)
    cop <- getAcop(family)
    V0 <- cop@V0(n, theta=theta)
    E <- matrix(rexp(n*asum), nrow=n, ncol=asum)
    X. <- E/V0 # psi^{-1}(U)
    ca <- c(0, cumsum(alpha))
    X <- sapply(2:length(ca), function(i)
                rowSums(X.[,(ca[i-1]+1):ca[i], drop=FALSE]))
    sapply(1:d, function(j) HbarL(X[,j], j=j, alpha=alpha, theta=theta,
                                  Rdist="Clayton"))
}
