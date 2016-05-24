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
########         BOOTSTRAP            FUNCTIONS         #########
#################################################################


# Find biases and standard errors for parameter estimates by parametric
# bootstrapping, and simple confidence intervals
# Simulate, many times, drawing samples from the estimated distribution, of
# the same size as the original data; re-estimate the parameters on the
# simulated data.  The distribution of the re-estimates around the estimated
# parameters is approximately the same as the distribution of the estimate
# around the true parameters.
# Invokes the estimating-equation MLE, but it would be easy to modify to
# use other methods.
# Confidence intervals (CI) are calculated for each parameter separately, using
# a simple pivotal interval (see, e.g., Wasserman, _All of Statistics_, section
# 8.3).  Confidence regions for combinations of parameters would be a tedious,
# but straightforward, extension.
# Inputs: distribution (as a list of the sort produced by tsal.fit), number of
#         bootstrap replicates, confidence level for confidence intervals,
#         distributional parameters (over-riding those of the distribution, if
#         one was given), fitting method (over-riding that used in the original
#         fit, if one was given), left-censoring threshold (over-riding, etc.)
# Outputs: structured list, containing the actual parameter settings used,
#          the estimated biases, the estimated standard errors, the lower
#          confidence limits, the upper confidence limits, the sample size, the
#          number of replicates, the confidence level, and the fitting method.
tsal.bootstrap.errors <- function(dist=NULL, reps=500, confidence=0.95,
    n=if(is.null(dist)) 1 else dist$n,
    shape=if(is.null(dist)) 1 else dist$shape,
    scale=if(is.null(dist)) 1 else dist$scale,
    q = if(is.null(dist)) tsal.q.from.shape(shape) else dist$q,
    kappa = if(is.null(dist)) tsal.kappa.from.ss(shape,scale) else dist$kappa,
    method = if(is.null(dist)) "mle.equation" else dist$method,
    xmin = if(is.null(dist)) 0 else dist$xmin) {
  # If we have both shape/scale and q/kappa parameters, the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  theta <- c(shape,scale,q,kappa)
  delta <- (1-confidence)/2 # probability of being in either tail, for CI
  CI.probs <- c(delta, 1-delta) # upper and lower probbilities for CI

  if(xmin == 0 || is.null(xmin))
  {
    Bootstrap <- replicate(reps, { x<-rtsal(n,shape,scale);
    est <- tsal.fit(x,xmin=0,method=method);
    c(est$shape, est$scale, est$q, est$kappa) })
  }else
  {
      Bootstrap <- replicate(reps, { x<-rtsal.tail(n,shape,scale,xmin=xmin);
                                 est <- tsal.fit(x,xmin=xmin,method=method);
                                 c(est$shape, est$scale, est$q, est$kappa) })
  }
  Bootstrap <- t(Bootstrap) # Make the columns the different variables
  Bootstrap.err <- t(apply(Bootstrap,1,"-",theta)) # Subtract theta from each
                                                   # row
  Bootstrap.bias <- apply(Bootstrap.err, 2, mean) # i.e. means of columns
  Bootstrap.se <- apply(Bootstrap.err, 2, sd) # std. devs. of columns
  Bootstrap.quantiles <- apply(Bootstrap.err,2,quantile,CI.probs) # quantiles
    # corresponding to upper and lower confidence limits, by columns
  CI.lower <- theta + Bootstrap.quantiles[1,] # lower confidence limits
  CI.upper <- theta + Bootstrap.quantiles[2,] # upper confidence limits
  namelist <- c("shape","scale","q","kappa")
  names(theta) <- namelist
  names(Bootstrap.bias) <- namelist
  names(Bootstrap.se) <- namelist
  names(CI.lower) <- namelist
  names(CI.upper) <- namelist
  return(list(originals=theta, bias=Bootstrap.bias, se=Bootstrap.se,
              confidence.interval.lower=CI.lower,
              confidence.interval.upper=CI.upper, sample.size=n,
              bootrap.replicates=reps, confidence=confidence, method=method,
              xmin=xmin))
}

# Estimate the total magnitude of a tail-sampled population
# Per Doug White's request
# Given that we have n samples from the tail of a distribution, i.e., only
# values >= xmin were retained, provide an estimate of the total magnitude
# (summed values) of the population
# Estimate the number of objects, observed and un-observed, as n/pr(X >= xmin)
# and then multiply by the mean
# Input: Distribution of the type produced by tsal.fit, distributional
#        parameters (over-riding the distribution if provided), number of
#        tail samples (over-riding the distribution if provided), xmin (over-
#        riding the one recorded in the distribution, if provided),
#        multiplier of size (if the base units of the data are not real units)
# Output: List, giving estimated total magnitude and estimated total population
#         size
tsal.total.magnitude <- function(dist=NULL, n=if(is.null(dist)) 1 else dist$n,
    shape=if(is.null(dist)) 1 else dist$shape,
    scale=if(is.null(dist)) 1 else dist$scale,
    q = if(is.null(dist)) tsal.q.from.shape(shape) else dist$q,
    kappa = if(is.null(dist)) tsal.kappa.from.ss(shape,scale) else dist$kappa,
    xmin = if(is.null(dist)) 0 else dist$xmin,
    mult = 1) {
  # If we have both shape/scale and q/kappa parameters, the latter over-ride.
  ss <- tsal.ss.from.qk(q,kappa)
  shape <- ss[1]
  scale <- ss[2]
  if(xmin == 0 || is.null(xmin))
  {
      p.tail <- ptsal(0,shape=shape,scale=scale,lower.tail=FALSE)
  }else
  {
      p.tail <- ptsal.tail(xmin,shape=shape,scale=scale,lower.tail=FALSE)
  }
  n.total <- n/p.tail
  mu <- tsal.mean(shape,scale)
  mag.total <- n.total * mu * mult
  totals <- list(magnitude.est = mag.total, count.est = n.total)
  return(totals)
}



