# From SamplerCompare, (c) 2010 Madeleine Thompson

# A slice sampler that samples along estimated eigenvectors of the
# underlying distribution.

univar.eigen.sample <- function(target.dist, x0, sample.size,
                                tuning=1, steps.out=100, cheat=FALSE) {
  ndim <- length(x0)
  X <- array(NA,c(sample.size,ndim))
  X[1,] <- x0
  nevals <- 0
  beta <- 0.05 # as with Roberts & Rosenthal
  burn.in <- 1
  nthin <- min(ndim, 25)
  
  # Eigendecomposition happens every decomp.freq univariate updates.
  # Never less than ten, and never such that it occurs more than a
  # hundred times in the simulation.

  decomp.freq <- max(ndim * floor(sample.size / 100), 10)

  # If we're cheating, obtain the covariance from the target distribution.

  if (cheat) {
    if (is.null(target.dist$cov))
      return(list(X=X[1,,drop=FALSE], evals=0))
    S.eig <- eigen(target.dist$cov)
  } else {
    S.eig <- NULL
  }

  # Sum of observations and scatter matrix.  Used to compute sample
  # covariance efficiently.

  obs.sum <- array(0, c(ndim, 1))
  obs.scatter <- array(0, c(ndim, ndim))
  scatter.N <- 0

  # Take sample.size*nthin univariate updates.

  for (obs in 2:(sample.size*nthin)) {

    # Every decomp.freq updates, update eigendecomposition if we're
    # still in the burn.in period (usually always) and we're not in
    # cheat mode.

    if (!cheat && (obs %% decomp.freq) == 0 && obs<nthin*sample.size*burn.in) {
      S.eig <- eigen(obs.scatter/scatter.N - tcrossprod(obs.sum/scatter.N))
    }

    # With prob. beta, pick a random direction.  Otherwise, pick a
    # random eigenvalue-eigenvector pair.

    if (runif(1) < beta || is.null(S.eig)) {
      v <- rnorm(ndim)
      v <- v / twonorm(v) * tuning
    } else {
      which.eig <- floor(1 + ndim * runif(1))
      v <- S.eig$vectors[,which.eig] * sqrt(abs(S.eig$values[which.eig]))
    }

    # Draw a slice level.

    y.slice <- target.dist$log.density(x0) - rexp(1)
    nevals <- nevals + 1

    # Position a unit interval [L,U] around zero, where coordinates
    # are relative to the current point in the direction of the
    # eigenvector with the square root of the eigenvalue being a unit
    # of distance.  Step out so that [L,U] is extends past the ends
    # of the slice.  See Neal '03, sec. 4 for discussion of algorithm.

    L <- -runif(1)
    U <- L + 1
    if (steps.out > 0) {
      L.y <- target.dist$log.density(x0+v*L) 
      U.y <- target.dist$log.density(x0+v*U) 
      nevals <- nevals + 2

      step <- 0
      while((L.y>y.slice || U.y>y.slice) && step < steps.out) {
        step <- step + 1
        if (runif(1)<0.5) {
          L <- L - 1
          L.y <- target.dist$log.density(x0+v*L) 
          nevals <- nevals + 1
        } else {
          U <- U + 1
          U.y <- target.dist$log.density(x0+v*U) 
          nevals <- nevals + 1
        }
      }
    }

    # Draw proposals from line along eigenvector, shrinking the
    # slice towards zero as prososals are rejected, until a proposal
    # is accepted.

    repeat {

      x1.offset <- runif(1, min=L, max=U)    # In scaled univariate coords.
      x1 <- x1.offset * v + x0               # In regular coords.

      y1 <- target.dist$log.density(x1)
      nevals <- nevals + 1

      if (y1>=y.slice)
        break

      if (x1.offset < 0)
        L <- x1.offset
      else
        U <- x1.offset
    }

    # Save every nthin updates into X and add the observation to the
    # sum and scatter matrix.

    if (obs %% nthin == 0) {
      X[obs/nthin,] <- x1
      obs.sum <- obs.sum + x1
      obs.scatter <- obs.scatter + tcrossprod(x1)
      scatter.N <- scatter.N + 1
    }
    x0 <- x1
  }
  return(list(X=X, evals=nevals, S.eig=S.eig))
}

attr(univar.eigen.sample, 'name') <- 'Univar Eigen'

# Curry out the cheat parameter so it can be passed as a function
# to compare.samplers.

cheat.univar.eigen.sample <- function(target.dist, x0, sample.size,
                                      tuning=1, steps.out=100) {
  univar.eigen.sample(target.dist, x0, sample.size, tuning, steps.out,
                      cheat=TRUE)
}

attr(cheat.univar.eigen.sample, 'name') <- 'Cheat Univar Eigen'
