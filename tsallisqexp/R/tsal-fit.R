# R functions for use with the q-exponential distribution of Tsallis

# Copyright (c) 2007, Cosma Shalizi, cshalizi@cmu.edu or cosma.shalizi@gmail.com

#####  This is free software; you can redistribute it and/or modify
#####  it under the terms of the GNU General Public License as published by
#####  the Free Software Foundation; either version 2 of the License, or
#####  (at your option) any later version.
#####
#####  This software is distributed in the hope that it will be useful,
#####  but WITHOUT ANY WARRANTY; without even the implied warranty of
#####  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#####  GNU General Public License for more details.
#####
#####  You should have received a copy of the GNU General Public License
#####  along with this software; if not, write to the Free Software
#####  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

# Please contact me at the above e-mail addresses for bug reports, questions,
# etc.  Please also contact me if you use this code in a scientific paper!

# This is a special case of the generalized Pareto distribution (type II), but
# arguably of independent interest.
# The distribution is defined through the upper cumulative or complementary
# distribution function, a.k.a. survival function,
### Pr(X>=x) = (1-(1-q)x/kappa)^(1/(1-q))
# It is convenient to introduce a re-parameterization
### shape = -1/(1-q)
### scale = shape*kappa
# which makes the relationship to the Pareto clearer, and eases estimation.
# Users can employ either parameterization, with shape/scale as the default,
# but q/kappa over-riding them.  (No warning is given if a user provides both
# sets of parameters.)
# For the derivation of the maximum likelihood estimators, please see the
# manuscript "Maximum Likelihood Estimation for q-Exponential (Tsallis)
# Distributions", which should accompany this file; if not it is available from
# http://bactra.org/research/tsallis-MLE/, along with the latest version of
# this code.  The manuscript is also available from
# http://arxiv.org/abs/math.ST/0701854

# Functions are divided into four parts:
# functions for the distribution itself, conforming to the R standards;
# functions to estimate the parameters;
# functions which illustrate the numerical accuracy of the implementation;
# functions to go from one parametrization to the other.


# Functions in this file:
# Distribution-related (per R standards):
# dtsal				Probability density
# ptsal				Cumulative probability
# qtsal				Quantiles
# rtsal				Random variate generation
# Parameter Estimation:
# tsal.loglik			Log-likelihood calculation
# tsal.fit			Estimate parameters; wrapper for detailed
#				methods.  Call this, not its subroutines
# tsal.mle.direct		Direct maximization of likelihood
# tsal.mle.equation		Find MLE by solving estimating equation; default
# tsal.est.shape.from.scale	MLE of shape parameter given scale parameter
# tsal.est.scale.from.shape	MLE of scale parameter given shape parameter
# tsal.curvefit			Find parameters by fitting a curve to the
#				empirical distribution function; avoid
# tsal.bootstrap.errors		Bootstrap standard errors, biases for MLE
# tsal.fisher			Fisher information matrix, for asymptotics
# tsal.mean			Calculate the expectation value
# tsal.total.magnitude		Total magnitude of a population (estimated)
# Implementation Testing:
# plot.tsal.quantile.transform	Illustrates relative numerical inaccuracy in
#				ptsal and qtsal, which should be inverses
# plot.tsal.LR.distribution	Calculates the log likelihood ratio for
#				estimated vs. fixed true parameters, and plots
#				it against the theoretical asymptotic
#				distribution (chi^2 with 2 d.f.).
# Censored Data:
# dtsal.tail			Probability density (tail-conditional)
# ptsal.tail			Cumulative probability (tail-conditional)
# qtsal.tail			Quantiles (tail-conditional)
# rtsal.tail			Random variate generation (from the tail)
# Parameter Conversion:
# tsal.shape.from.q		Get shape parameter from q
# tsal.scale.from.qk		Get scale parameter from q, kappa
# tsal.q.from.shape		Get q from from shape parameter
# tsal.kappa.from.ss		Get kappa from shape, scale
# tsal.ss.from.qk		Get shape, scale from q, kappa (as pairs)
# tsal.qk.from.ss		Get q, kappa from shape, scale (as pairs




#################################################################
########         PARAMETER ESTIMATION FUNCTIONS         #########
#################################################################


# Calculate log-likelihood
# If a left-censoring threshold is set, and any of the data are below it,
# this should produce NA, probably with an error message
# Input: Data vector, distributional parameters, left-censoring threshold
# Output: log likelihood
tsal.loglik <- function(x, shape, scale, q=tsal.q.from.shape(shape),
                         kappa=tsal.kappa.from.ss(shape,scale),
                         xmin=0) {
  # If we have both shape/scale and q/kappa parameters, the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  if(xmin == 0 || is.null(xmin))
  {
      l <- sum(dtsal(x,shape=shape,scale=scale,log=TRUE))
  }else
  {
      l <- sum(dtsal.tail(x,shape=shape,scale=scale,xmin=xmin,log=TRUE))
  }
  return(l)
}



# Fit a q-exponential to data
# Wrapper for the actual methods:
#     tsal.fit.mle.equation (solve maximum likelihood estimating equations);
#     tsal.fit.mle.direct (numerical likelihood maximization); and
#     tsal.fit.leastsquares (least-squares curve-fitting to the empirical
#                            distribution);
# prettying up the results in all cases.
# Users should call this function with the optional method argument, NOT the
# subsidiary functions.  Allowed values are the extensions of the functions
# listed above, plus "curvefit" as an alias to leastsquares.
# The likelihood methods are preferred, since the MLE is consistent and
# asymptotically unbiased and efficient.  The default is to solve the
# transcendental estimating equations numerically, as it is by far the fastest
# method.  Occasionally, however, R's solver gives weird results, with huge
# shape/scale parameters, and then direct optimization of the likelihood can be
# superior.
# The curve-fitting method is disrecommended, and included solely for purposes
# of "backwards compatibility" with the practice of the statistical physics
# community.
# Input: Data, left-censoring threshold, method tag, optional arguments to
#        optimization routines (if invoked)
# Output: A list having as components the shape/scale and q/kappa parameters,
#         plus some information about the fit and fitting controls
tsal.fit <- function(x, xmin=0, method=c("mle.equation", "mle.direct", "leastsquares"),...) {
    
    method <- match.arg(method, c("mle.equation", "mle.direct", "leastsquares"))
    
  switch(method,
    mle.equation = { scale.est <- tsal.mle.equation(x,xmin) },
    mle.direct = { scale.est <- tsal.mle.direct(x,xmin,...) },
    curvefit = { scale.est <- tsal.fit.leastsquares(x,...) },
    leastsquares = { scale.est <- tsal.fit.leastsquares(x,...) },
    { cat("Unknown method ", method, "\n");  scale.est <- NA }
  )
  if (!is.na(scale.est)) {
    shape.est <- tsal.est.shape.from.scale(x,scale.est,xmin=xmin)
    qk <- tsal.qk.from.ss(shape.est,scale.est)
    q.est <- qk[1]
    kappa.est <- qk[2]
    loglik <- tsal.loglik(x,shape=shape.est,scale=scale.est,xmin=xmin)
    fit <- list(type="tsal", q=q.est, kappa=kappa.est, shape=shape.est,
                scale=scale.est, loglik = loglik, n = length(x), xmin=xmin,
                method=method)
  } else {
     fit <- NA
  }
  return(fit)   
}

# Fit a q-exponential to data by solving the likelihood estimating equations
# Returns ONLY the scale parameter, users should always call tsal.mle, which
# does everything else
# Uses the basic R routine for finding the roots of a one-variable equation,
# which requires an initial bracket on the values; if the defaults don't work,
# I have not been clever enough to come up with a good bracketing procedure,
# and several crude ones give very bad results, so I give up and pass the
# problem to the numerical optimnization routine tsal.mle.direct, which is
# slow but always seems to give reasonable-looking results.
# Input: Data vector, left-censoring threshold
# Output: estimate of scale parameter
tsal.mle.equation <- function(x,xmin=0) {
  n <- length(x)
  # Set up the function so that it handles the case with xmin possibly != 0
  # This form derived from the equation in the paper by multiplying
  # numerators and denominators by the scale factor, and multiplying
  # through by n, all in the name of minimizing division
  g <- function(scale) {
    shape <- tsal.est.shape.from.scale(x,scale,xmin)
    z <- x+scale
    y <- x/z # should work via recycling rule
    return(n - (shape+1)*sum(y) + n*shape*xmin/(scale+xmin))
  }
  # Initial guesses as to brackets of the root
  # For positive shape, or q > 1
  lower.p <- 1e-12 # Zero doesn't quite work here
  upper.p <- max(x)
  # For negative shape, or q < 1
  lower.m <- -100*max(x)
  upper.m <- -max(x)-1 # Has to be strictly less than -max(x)
  # Check whether the brackets make sense
  # Don't evaluate the function multiple times
  g.low.p <- g(lower.p)
  g.up.p <- g(upper.p)
  g.low.m <- g(lower.m)
  g.up.m <- g(upper.m)
  tried.p <- FALSE
  tried.m <- FALSE
  # The brackets make sense, solve
  if ((sign(g.low.p) != sign(g.up.p)) && (g.up.p!=0) && (g.low.p!=0)) {
    scale.est.p <- uniroot(g,lower=lower.p,upper=upper.p)$root
    tried.p <- TRUE
  }
  if ((sign(g.low.m) != sign(g.up.m)) && (g.up.m!=0) && (g.low.m!=0)) {
    scale.est.m <- uniroot(g,lower=lower.m,upper=upper.m)$root
    tried.m <- TRUE
  }
  # If the brackets don't work, don't bother trying to fix them, just optimize
  # the likelihood
  if ((!tried.p) && (!tried.m)) {
    return(tsal.mle.direct(x,xmin))
  }
  # Otherwise, if at least one bracket worked, pick the best
  if ((tried.p) && (!tried.m)) { scale.est <- scale.est.p }
  if ((!tried.p) && (tried.m)) { scale.est <- scale.est.m }
  if ((tried.p) && (tried.m)) {
    # We don't ever seem to go here, but I'm not sure we can't
    shape.est.p <- tsal.est.shape.from.scale(x,scale.est.p,xmin)
    shape.est.m <- tsal.est.shape.from.scale(x,scale.est.m,xmin)
    lp <- tsal.loglik(x,shape.est.p,scale.est.p,xmin=xmin)
    lm <- tsal.loglik(x,shape.est.m,scale.est.m,xmin=xmin)
    if (lm > lp) { scale.est <- scale.est.m } else { scale.est <- scale.est.p }
  }
  # Otherwise, use those brackets to find a root
  return(scale.est)
}

# Fit a q-exponential to data by numerically maximizing the likelihood
# Returns ONLY the scale parameter; users should always call tsal.mle, which
# does everything else
# Input: data, left-censoring threshold, optional arguments to maximizer
# Output: estimate of scale parameter
tsal.mle.direct <- function(x,xmin,...) {
  l <- function(t) {-tsal.loglik(x,shape=t[1],scale=t[2],xmin=xmin)}
  # Starting parameter values for positive shape, or q>1
  t0.plus <- c(1,max(x)+1)
  # Starting parameter values for negative shape, or q<1
  t0.minus <- c(-1,-max(x)-1)
  est.plus <- optim(par=t0.plus, fn=l,gr=NULL,...)
  est.minus <- optim(par=t0.minus, fn=l,gr=NULL,...)
  scale.est <- est.plus$par[2]
  if (est.minus$value < est.plus$value) { scale.est <- est.minus$par[2] }
  return(scale.est)
}

# Compute MLE of a q-exponential's shape, assuming scale is known
# Input: Data vector, scale parameter, left-censoring threshold
# Output: Real-valued shape estimate
tsal.est.shape.from.scale <- function(x,scale,xmin=0) {
  n <- length(x)
  z <- plus((scale+x)/(scale+xmin))
  shape.est <- n/sum(log(z))
  return(shape.est)
}

# Compute MLE of a q-exponential's scale, assuming shape is known
# Input: Data vector, scale, left-censoring threshold
# Output: Real-valued scale estimate
tsal.est.scale.from.shape <- function(x,shape,xmin=0) {
  n <- length(x)
  # Set up the function so that it handles the case with xmin possible != 0
  # This form derived from the equation in the paper by multiplying
  # numerators and denominators by the scale factor, and multiplying
  # through by n, all in the name of minimizing division
  f <- function(scale) { 
    z <- scale+x
    y <- x/z # should work via recycling rule
    return(n-(shape+1)*sum(y) + n*shape*xmin/(scale+xmin))
  }
  lower <- 0
  upper <- max(x)*(shape+1)
  scale.est <- uniroot(f,lower,upper)$root
  return(scale.est)
}

# Estimate the parameters of a q-exponential distribution by least-squares
# curve-fitting to the empirical cumulative distribution function
# Returns ONLY the scale parameter, users should always call tsal.mle, which
# does everything else
# THEORY:
# The survival function is
### Pr(X>=x) = (1+x/scale)^(-shape)
# therefore
### log(Pr(X>=x)) + shape*log(1+x/scale) = 0
# Writing S_n for the empirical survival function, the least-squares solution
# to
### sum_{i}{(log(S_n(x_i)) + shape*log(1+x_i/scale))^2}
# COULD be used to estimate the parameters.
# This seems to be the current (2007) best practice in statistical physics.
# However, it is, or ought to be, well-known that the analogous estimator for
# the Pareto distribution is highly inaccurate, and inferior to the maximum
# likelihood estimator.  The same is true for q-exponential distributions.
# This function is only for comparison and "backward compatibility", and I
# strongly recommend against using it.
# The default optimization algorithm is Nelder-Mead.  The search algorithm and
# other aspects of the optimization can be changed by passing in optional
# arguments.  (See the R help on "optim" for details.)
# Inputs: Data vector, optional arguments for the optimizer
# Outputs: Real-valued estimate of the scale parameter
tsal.fit.leastsquares <- function(x, ...) {
  n <- length(x)
  Fn <- ecdf(x)
  k <- knots(Fn)
  Sn <- 1 + (1/n) - Fn(k)  # Avoids the zero value at the last point
  f <- function(theta) { sum((log(Sn) + theta[1]*log(1+k/theta[2]))^2) }
  theta_0 <- c(1,1) # Arbitrary, but we've got to start somewhere!
  scale.est <- optim(par=theta_0,fn=f,gr=NULL,...)$par[2]
  return(scale.est)
}


# Calculate the Fisher information matrix, for asymptotic variances and
# covariances of the maximum likelihood estimates of shape and scale
# First row/column corresponds to shape, second to scale
# Convergence to the asymptotic normal distribution can be slow, so for limited
# data you should bootstrap
# Note that this function ONLY works with the shape-scale parameterization
# Inputs: shape, scale, left-censoring threshold
# Outputs: a matrix
tsal.fisher <- function(shape, scale, xmin=0) {
    # If xmin > 0, then effectively the scale = scale+xmin (see paper)
    scale <- scale+xmin
    shape.shape <-  scale^(-2)
    shape.scale <- -1/(scale*(shape+1))
    scale.scale <- shape/(scale^2 * (shape+2))
    fisher <- rbind(c(shape.shape,shape.scale),c(shape.scale,scale.scale))
    return(as.matrix(fisher))
}


