
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:                     DESCRIPTION:
# .cov.shrink.tawny
# .getCorFilter.Shrinkage
# .cov.sample.tawny
# .cov.prior.cc
# .cov.prior.identity
# .cor.mean.tawny
# .shrinkage.intensity
# .shrinkage.p
# .shrinkage.r
# .shrinkage.c
################################################################################


# Rmetrics: 
#   Note that tawny is not available on Debian as of 2009-04-28. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.
#   We also made modifications for tailored usage with Rmetrics. 


# Package: tawny
# Title: Provides various portfolio optimization strategies including
#    random matrix theory and shrinkage estimators
# Version: 1.0
# Date: 2009-03-02
# Author: Brian Lee Yung Rowe
# Maintainer: Brian Lee Yung Rowe <tawny-help@muxspace.com>
# License: GPL-2
#
# Perform shrinkage on a sample covariance towards a biased covariance
#
# This performs a covariance shrinkage estimation as specified in Ledoit 
# and Wolf. Using within the larger framework only requires using the 
# getCorFilter.Shrinkage function, which handles the work of constructing 
# a shrinkage estimate of the covariance matrix of returns (and consequently 
# its corresponding correlation matrix).


# ------------------------------------------------------------------------------


.cov.shrink.tawny <- 
function(returns, sample = NULL, prior.fun = .cov.prior.cc, ...)
{
    # Shrink the sample covariance matrix towards the model covariance 
    #   matrix for the given time window.
    # model - The covariance matrix specified by the model, e.g. single-index, 
    #   Barra, or something else
    # sample - The sample covariance matrix. If the sample covariance is null, 
    #   then it will be computed from the returns matrix
    
    # Example
    #   S.hat <- .cov.shrink.tawny(ys)

    # if (.loglevel.tawny() > 0) cat("Shrinking covariance for",last(index(returns)),"\n")
    if (is.null(sample)) { S <- .cov.sample.tawny(returns) }
    else { S <- sample }
    
    T <- nrow(returns)
    # F <- .cov.prior.cc(S)
    F <- prior.fun(S, ...)
    k <- .shrinkage.intensity(returns, F, S)
    d <- max(0, min(k/T, 1))
    
    if (.loglevel.tawny() > 0) cat("Got intensity k =", k,
        "and coefficient d =",d,"\n")
    
    S.hat <- d * F + (1 - d) * S
    S.hat
}

# ------------------------------------------------------------------------------


.getCorFilter.Shrinkage <- 
function(prior.fun = .cov.prior.cc, ...)
{
    # Return a correlation matrix generator that is compatible with the 
    # portfolio optimizer
    
    # Example
    #   ws <- optimizePortfolio(ys, 100, .getCorFilter.Shrinkage())
    #   plotPerformance(ys,ws)

function(h) return(cov2cor(.cov.shrink.tawny(h, prior.fun=prior.fun, ...)))
}


# ------------------------------------------------------------------------------


.cov.sample.tawny <- 
function(returns)
{
    # Calculate the sample covariance matrix from a returns matrix
    # Returns a T x N returns matrix 
    # p.cov <- .cov.sample.tawny(p)

    # X is N x T
    T <- nrow(returns)
    X <- t(returns)
    ones <- rep(1,T)
    S <- (1/T) * X %*% (diag(T) - 1/T * (ones %o% ones) ) %*% t(X)
    S
}

# ------------------------------------------------------------------------------


.cov.prior.cc <- 
function(S)
{
    # Constant correlation target
    # S is sample covariance

    r.bar <- .cor.mean.tawny(S)
    vars <- diag(S) %o% diag(S)
    F <- r.bar * (vars)^0.5
    diag(F) <- diag(S)
    return(F)
}


# ------------------------------------------------------------------------------


.cov.prior.identity <- 
function(S)
{
    # This returns a covariance matrix based on the identity (i.e. no 
    # correlation)
    
    # S is sample covariance

    return(diag(nrow(S)))
}


# ------------------------------------------------------------------------------


.cor.mean.tawny <- 
function(S)
{
    # Get mean of correlations from covariance matrix
    
    N <- ncol(S)
    cors <- cov2cor(S)
    2 * sum(cors[lower.tri(cors)], na.rm=TRUE) / (N^2 - N)
}


# ------------------------------------------------------------------------------
  
                                                                             
.shrinkage.intensity <- 
function(returns, prior, sample)
{
    # Calculate the optimal shrinkage intensity constant
    # returns : asset returns T x N
    # prior : biased estimator

    p <- .shrinkage.p(returns, sample)
    
    r <- .shrinkage.r(returns, sample, p)
    c <- .shrinkage.c(prior, sample)
    (p$sum - r) / c
}


# ------------------------------------------------------------------------------


.shrinkage.p <- 
function(returns, sample)
{
    # Sum of the asymptotic variances
    # returns : T x N - Matrix of asset returns
    # sample : N x N - Sample covariance matrix
    # Used internally.
    # S <- .cov.sample.tawny(ys)
    # ys.p <- .shrinkage.p(ys, S)

    T <- nrow(returns)
    N <- ncol(returns)
    ones <- rep(1,T)
    means <- t(returns) %*% ones / T
    z <- returns - matrix(rep(t(means), T), ncol=N, byrow=TRUE)
    
    term.1 <- t(z^2) %*% z^2
    term.2 <- 2 * sample * (t(z) %*% z)
    term.3 <- sample^2
    phi.mat <- (term.1 - term.2 + term.3) / T
    
    phi <- list()
    phi$sum <- sum(phi.mat)
    phi$diags <- diag(phi.mat)
    phi
}


# ------------------------------------------------------------------------------


.shrinkage.r <- 
function(returns, sample, pi.est)
{
    # Estimation for rho when using a constant correlation target
    # returns : stock returns
    # market : market returns 
    # Example
    #   S <- .cov.sample.tawny(ys)
    #   ys.p <- .shrinkage.p(ys, S)
    #   ys.r <- .shrinkage.r(ys, S, ys.p)

    N <- ncol(returns)
    T <- nrow(returns)
    ones <- rep(1,T)
    means <- t(returns) %*% ones / T
    z <- returns - matrix(rep(t(means), T), ncol=N, byrow=TRUE)
    r.bar <- .cor.mean.tawny(sample)
    
    # Asymptotic covariance estimator
    term.1 <- t(z^3) %*% z
    term.2 <- diag(sample) * (t(z) %*% z)
    term.3 <- sample * (t(z^2) %*% matrix(rep(1,N*T), ncol=N))
    # This can be simplified to diag(sample) * sample, but this expansion is
    # a bit more explicit in the intent (unless you're an R guru)
    term.4 <- (diag(sample) %o% rep(1,N)) * sample
    script.is <- (term.1 - term.2 - term.3 + term.4) / T
    
    # Create matrix of quotients
    ratios <- (diag(sample) %o% diag(sample)^-1)^0.5
    
    # Sum results
    rhos <- 0.5 * r.bar * (ratios * script.is + t(ratios) * t(script.is))
    
    # Add in sum of diagonals of pi
    sum(pi.est$diags, na.rm = TRUE) 
    + sum(rhos[lower.tri(rhos)], na.rm = TRUE) 
    + sum(rhos[upper.tri(rhos)], na.rm = TRUE)
}


# ------------------------------------------------------------------------------

.shrinkage.c <- 
function(prior, sample)
{
    # Misspecification of the model covariance matrix

    squares <- (prior - sample)^2
    sum(squares, na.rm = TRUE)
}


# ------------------------------------------------------------------------------


.loglevel.tawny <- 
function (new.level = NULL) 
{
    if (!is.null(new.level)) {
        options(log.level = new.level)
    }
    
    if (is.null(getOption("log.level"))) {
        return(0)
    }
    
    return(getOption("log.level"))
}



################################################################################

