# From SamplerCompare, (c) 2010 Madeleine Thompson

# metropolis.R contains implementations of three variations on the
# Metropolis sampler: Metropolis with spherical proposals, Metropolis
# with univariate proposals, and adaptive Metropolis.

# multivariate.metropolis.sample simulates a target distribution
# with multivariate spherical Metropolis.  See
# ?multivariate.metropolis.sample for more information.

multivariate.metropolis.sample <-
    function(target.dist, x0, sample.size, tuning=1) {
  stopifnot(target.dist$ndim==length(x0))
  X <- array(NA,c(sample.size,target.dist$ndim))
  nevals <- 1
  y0 <- target.dist$log.density(x0)     # log density at current state
  rej <- 0                              # number of rejected proposals

  for (obs in 1:(sample.size*target.dist$ndim)) {
    x1 <- rnorm(target.dist$ndim, x0, tuning)   # proposal
    nevals <- nevals + 1
    y1 <- target.dist$log.density(x1)           # log density at proposal
    if (runif(1) < exp(y1-y0)) {                # accept proposal?
      x0 <- x1
      y0 <- y1
    } else {
      rej <- rej + 1
    }
    if ((obs %% target.dist$ndim) == 0)         # keep every ndim observations
      X[obs/target.dist$ndim,] <- x0
  }

  reject.rate <- rej/(sample.size*target.dist$ndim)
  return(list(X=X, reject.rate=reject.rate, evals=nevals))
}

attr(multivariate.metropolis.sample, 'name') <- 'Multivariate Metropolis'

# univar.metropolis.sample simulates a target distribution with
# Metropolis with single coordinate updates.  See ?univar.metropolis.sample
# for more information.

univar.metropolis.sample <- function(target.dist, x0, sample.size, tuning=1) {
  stopifnot(target.dist$ndim==length(x0))
  X <- array(NA,c(sample.size,target.dist$ndim))
  X[1,] <- x0
  nevals <- 1
  y0 <- target.dist$log.density(x0)         # log density at current state
  rej <- 0                                  # number of rejected proposals

  for (obs in 2:sample.size) {              # iterate through sample.size
    x0 <- X[obs-1,,drop=TRUE]

    for (i in 1:target.dist$ndim) {         # iterate through coordinates
      x1 <- x0
      x1[i] <- rnorm(1, x0[i], tuning)      # propose change in one coord.
      nevals <- nevals + 1
      y1 <- target.dist$log.density(x1)     # log-density at proposal
      if (runif(1) < exp(y1-y0)) {          # accept proposal?
        x0 <- x1
        y0 <- y1
      } else {
        rej <- rej + 1
      }
    }
    X[obs,] <- x0                           # coordinates updated, save state
  }

  reject.rate <- rej/(sample.size*target.dist$ndim)
  return(list(X=X, reject.rate=reject.rate, evals=nevals))
}

attr(univar.metropolis.sample, 'name') <- 'Univariate Metropolis'


# adaptive.metropolis.sample simulates a target distribution with
# the adaptive Metropolis algorithm of Roberts and Rosenthal (2009).
# More information can be found in ?adaptive.metropolis.sample.

adaptive.metropolis.sample <-
    function(target.dist, x0, sample.size, tuning=0.1, beta=0.05, burn.in=0.2) {
  stopifnot(target.dist$ndim==length(x0))
  ndim <- target.dist$ndim
  X <- array(NA, c(sample.size, ndim))
  nevals <- 1
  y0 <- target.dist$log.density(x0)     # log density at current state
  rej <- 0                              # number of rejected proposals

  # Sum of observations and scatter matrix of observations with
  # respect to zero.  Used to efficiently estimate covariance.

  obs.sum <- array(0, c(ndim, 1))
  obs.scatter <- array(0, c(ndim, ndim))

  # Most recently computed sample covariance.
  sample.cov <- NULL

  # Cholesky factor of empirical component of proposal variance
  prop.R <- NULL

  for (obs in 1:(sample.size*ndim)) {

    # Draw a proposal from a mixture of Gaussians, one spherical
    # and one proportional to the empirical distribution.  See eqn.
    # 2.1 of Roberts & Rosenthal (2009).

    if (is.null(prop.R) || runif(1) < beta) {
      x1 <- rnorm(ndim, x0, tuning/sqrt(ndim))
    } else {
      x1 <- t(prop.R) %*% rnorm(ndim) + x0
    }

    # Compute log density and test for proposal acceptance.

    nevals <- nevals + 1
    y1 <- target.dist$log.density(x1)           # log density at proposal
    if (runif(1) < exp(y1-y0)) {                # accept proposal?
      x0 <- x1
      y0 <- y1
    } else {
      rej <- rej + 1
    }

    # Update scatter and sample sum.

    obs.sum <- obs.sum + x0
    obs.scatter <- obs.scatter + tcrossprod(x0)

    # Every ndim observations, keep one.

    if ((obs %% ndim) == 0) {
      X[obs/ndim,] <- x0
    }

    # Every 10*ndim observations, update proposal covariance.  This
    # is O(ndim^3), and drawing a proposal is O(ndim^2), so this rate
    # hopefully keeps proposal update computation requirements reasonable.

    if ((obs %% (10*ndim)) == 0 && obs<sample.size*burn.in) {
      sample.cov <- obs.scatter/obs - tcrossprod(obs.sum/obs)
      prop.R <- try(2.38/sqrt(ndim) * chol(sample.cov), silent=TRUE)
      if (!is.matrix(prop.R))  # the sample has too few unique values
        prop.R <- NULL
    }
  }

  reject.rate <- rej/(sample.size*target.dist$ndim)
  return(list(X=X, reject.rate=reject.rate, evals=nevals,
              sample.cov=sample.cov))
}

attr(adaptive.metropolis.sample, 'name') <- 'Adaptive Metropolis'
