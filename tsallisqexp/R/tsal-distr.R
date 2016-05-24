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
########       DISTRIBUTION-RELATED FUNCTIONS           #########
#################################################################


# Calculate the probability density
# Input: vector of data values, distributional parameters, left-censoring
#        threshold, log flag
# Output: vector of (log) densities
dtsal <- function(x, shape=1,scale=1, q=tsal.q.from.shape(shape),
                  kappa=tsal.kappa.from.ss(shape,scale),
                  log=FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  # Under censoring, pass off to the tail version
  z <- plus(1+x/scale)
  if (log)
  {
      d <- -(shape+1)*log(z) + log(shape/scale)
      d[is.infinite(d)] <- -Inf
  }else
  {
      d <- (shape/scale)*(z^(-shape-1))
      d[is.infinite(d)] <- 0
  }
  return(d)
}

# Calculate the cumulative distribution function
# Input: vector of data values, distributional parameters, left-censoring
#        threshold, usual flags
# Output: vector of (log) probabilities
ptsal <- function(x, shape=1, scale=1, q=tsal.q.from.shape(shape),
                  kappa=tsal.kappa.from.ss(shape,scale),
                  lower.tail=TRUE, log.p=FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]

  z <- plus(1+x/scale)
  if ((log.p) && (!lower.tail)) { p <- -shape*log(z) }
  if ((log.p) && (lower.tail)) { p <- log(1-(z^(-shape))) }
  if ((!log.p) && (!lower.tail))
  {
      p <- z^(-shape)
      p[p > 1] <- 1
  }
  if ((!log.p) && (lower.tail))
  {
      p <- 1 - z^(-shape)
      p[p < 0] <- 0
  }
  return(p)
}

# Calculate quantiles
# Input: vector of p-values, distributional parameters, left-censoring
#        threshold, usual flags
# Output: vector of quantile locations
qtsal <- function(p,  shape=1, scale=1, q=tsal.q.from.shape(shape),
                  kappa=tsal.kappa.from.ss(shape,scale),
                  lower.tail=TRUE, log.p=FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]

  if (log.p) { p <- exp(p) }
  if (lower.tail) { p <- 1-p }
  # The upper quantile function is given by
  # (scale)(p^{-1/shape} - 1) = x
  quantiles <- scale*(-1 + (p^(-1/shape)))
  return(quantiles)
}

# Generate random variates
# Input: integer length, distributional parameters, left-censoring threshold
# Output: vector of reals
rtsal <- function(n, shape=1, scale=1, q=tsal.q.from.shape(shape),
                 kappa=tsal.kappa.from.ss(shape,scale)) {
  # If we have both shape/scale and q/kappa parameters, the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]

  # Apply the transformation method
  ru <- runif(n)
  r <- qtsal(ru,shape,scale)
  return(r)
}



# Calculate the expectation value
# Input: Distributional parameters
# Output: The mean of that distribution
tsal.mean <- function(shape, scale, q=tsal.q.from.shape(shape),
kappa=tsal.kappa.from.ss(shape,scale)) {
    # If we have both shape/scale and q/kappa parameters, the latter over-ride.
    ss <- tsal.ss.from.qk(q,kappa)
    shape <- ss[1]
    scale <- ss[2]
    mu <- scale/(shape-1)
    return(mu)
}


