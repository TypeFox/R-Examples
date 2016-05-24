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
# tsal.loglik		Log-likelihood calculation
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
########          CENSORED DATA FUNCTIONS               #########
#################################################################

# Users will have little reason to call these directly, instead of using
# the higher-level ones with the optional xmin argument

# Calculate the probability density, conditional on being in the right tail
# Users should call dtsal with the xmin argument, not this
# Outputs NA at values below the threshold
# Input: vector of data values, distributional parameters, log flag
# Output: vector of (log) densities
dtsal.tail <- function(x, shape=1,scale=1, q=tsal.q.from.shape(shape),
                       kappa=tsal.kappa.from.ss(shape,scale), xmin=0,
                       log=FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  # Default to the non-tail version when xmin==0
  # but this function should never be called with xmin==0
  if(xmin==0) { return(dtsal(x, shape, scale, q, kappa, log)) }
  z <- plus(1+x/scale)
  C <- 1/ptsal(xmin,shape,scale,lower.tail=FALSE)
  C.log <- -ptsal(xmin,shape,scale,lower.tail=FALSE,log.p=TRUE)
  
  # Possible later modification: Check whether C is so small that ONLY C.log
  # should be reliably used.  ("slide-rule" subroutine?)
  if (log)
  {
      f <- function(z) { C.log - (shape+1)*log(z) + log(shape/scale) }
      d <- ifelse(x < xmin, -Inf, f(z))
  }else
  {
      f <- function(z) {C*(shape/scale)*(z^(-shape-1)) }
      d <- ifelse(x < xmin, 0, f(z))
  }
  return(d)
}

# Calculate the cumulative distribution function, conditional on being in the
# right tail
# Users should call ptsal with the xmin argument, not this
# Outputs NA at values below the threshold
# xmin argument
# Input: vector of data values, distributional parameters, usual flags
# Output: vector of (log) probabilities
ptsal.tail <- function(x, shape=1, scale=1, q=tsal.q.from.shape(shape),
                  kappa=tsal.kappa.from.ss(shape,scale), xmin=0,
                  lower.tail=TRUE, log.p=FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  # Default to the non-tail version when xmin==0
  # but this function should never be called with xmin==0
  if(xmin==0) { return(ptsal(x, shape, scale, q, kappa, lower.tail, log.p)) }
  C <- 1/ptsal(xmin,shape,scale,lower.tail=FALSE)
  C.log <- -ptsal(xmin, shape, scale, lower.tail=FALSE, log.p=TRUE)
  z <- plus(1+x/scale)
  
  if ((log.p) && (!lower.tail))
  {
      f <- function(z) { C.log -shape*log(z) }
      p <- ifelse(x < xmin, 0, f(z))
  }
  if ((log.p) && (lower.tail))
  {
      f <- function(z) { log(1-(C*z^(-shape))) }
      p <- ifelse(x <= xmin, -Inf, f(z))
  }
  if ((!log.p) && (!lower.tail))
  {
      f <- function(z) { C*z^(-shape) }
      p <- ifelse(x < xmin, 1, f(z))
  }
  if ((!log.p) && (lower.tail))
  {
      f <- function(z) { 1 - C*z^(-shape) }
      p <- ifelse(x < xmin, 0, f(z))
  }
  
  return(p)
}

# Calculate quantiles, conditional on being in the right tail
# Users should call qtsal with the xmin argument, not this
# The tail-conditional quantile function is related to the usual one by
# Q_tail(p) = Q(1 - p S(xmin))
# where Q is the usual quantile function and S is the UNCONDITIONAL
# survival function
# Input: vector of p-values, distributional parameters, threshold, usual flags
# Output: vector of quantile locations
qtsal.tail <- function(p,  shape=1, scale=1, q=tsal.q.from.shape(shape),
                  kappa=tsal.kappa.from.ss(shape,scale), xmin=0,
                  lower.tail=TRUE, log.p=FALSE) {
  # If we have both shape/scale and q/kappa parameters,
  # the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  # Default to the non-tail version when xmin=0
  # but this function should never be called with xmin==0
  if(xmin==0) { return(qtsal(p, shape, scale, q, kappa, lower.tail, log.p)) }
  # Transform the probabilities, if necessary
  if (log.p) { p <- exp(p) }
  if (lower.tail) { p <- 1-p }
  
  # The right-tail-conditional quantile function, Q_tail, is related to the
  # usual one via
  # Q_tail(p) = Q(1 - p S(xmin))
  # where S is the UNCONDITIONAL survival function
  # Let's compute that constant; only let's do it's log, because it can be very
  # small
  C.log <- ptsal(xmin, shape, scale, lower.tail=FALSE, log.p=TRUE)
  # Now multiply it by the probabilities
  p.multiplied <- exp(C.log + log(p))
  # This possibly undoes some of the work under log.p above
  # Now invoke the unconditional quantile function
  quantiles <- qtsal(1-p.multiplied,shape,scale)
  return(quantiles)
}

# Generate random variates from the right tail
# Users should call rtsal with the xmin argument, not this
# Input: integer length, distributional parameters, threshold
# Output: vector of reals
rtsal.tail <- function(n, shape=1, scale=1, q=tsal.q.from.shape(shape),
                       kappa=tsal.kappa.from.ss(shape,scale), xmin=0) {
  # If we have both shape/scale and q/kappa parameters, the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  # Default to the non-tail version when xmin==0
  # but this function should never be called with xmin==0
  if(xmin==0) { return(rtsal(n, shape, scale)) }
  # Apply the transformation method
  ru <- runif(n)
  r <- qtsal.tail(ru,shape,scale,xmin=xmin)
  return(r)
}
