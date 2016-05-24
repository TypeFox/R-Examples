# Two variations of the hyperrectangle method that adapt the
# orientation of the slice based on the eigendecomposition of the
# covariance.
#
# oblique.hyperrect.sample: estimates the covariance with the MLE
#
# cheat.oblique.hyperrect.sample: cheats and uses target.dist$cov

oblique.hyperrect.sample <- function(target.dist, x0, sample.size, tuning=1,
                                     edge.scale=5, cheat=FALSE) {
  ndim <- length(x0)
  X <- array(NA,c(sample.size,ndim))
  X[1,] <- x0
  nevals <- 0
  beta <- 0.05 # as with Roberts & Rosenthal
  burn.in <- 1
  
  # Eigendecomposition happens every decomp.freq univariate updates.
  # Never less than ten, and never such that it occurs more than a
  # hundred times in the simulation.

  decomp.freq <- max(floor(sample.size / 100), 10)

  # If we're cheating, obtain the covariance from the target distribution.
  
  if (cheat) {
    if (is.null(target.dist$cov))
      return(list(X=X[1,,drop=FALSE], evals=0))
    S.eig <- eigen(target.dist$cov)
  } else {
    S.eig <- NULL
  }

  for (obs in 2:sample.size) {

    # Every decomp.freq updates, update eigendecomposition if we're
    # still in the burn.in period (usually always) and we're not in
    # cheat mode.

    if (!cheat && (obs %% decomp.freq) == 0 && obs<sample.size*burn.in) {
      S.eig <- eigen(cov(X[1:(obs-1),,drop=FALSE]))
    }

    # In the coordinate system defined by the eigenvectors and
    # eigenvalues with origin at the current point, position a unit
    # hypercube so that the origin is uniformly distributed with
    # respect to its corners.

    L <- -1 * runif(ndim)
    U <- L + 1

    # With probability beta, approximate the slice with a hypercube
    # with edge length [tuning].  Otherwise, approximate it with the
    # eigenvalues and eigenvectors.

    if (runif(1)<beta || is.null(S.eig)) {
      vals <- rep(tuning, ndim)
      vecs <- diag(1, nrow=ndim)
    } else {
      vals <- S.eig$values
      vecs <- S.eig$vectors
    }

    # Draw a slice level.

    y.slice <- target.dist$log.density(x0) - rexp(1)
    nevals <- nevals + 1

    # Draw proposals till one is accepted.

    repeat {
      # Draw a proposal from the approximation.  wt is coordinates
      # in the transformed space.  v is an offset from the current
      # point in the original space.  x1 is a proposal in the original space.

      wt <- runif(ndim, min=L, max=U)
      v <- as.numeric(vecs %*% (edge.scale * wt * vals))
      x1 <- x0 + v

      # Accept the proposal if it is in the slice.

      y1 <- target.dist$log.density(x1)
      nevals <- nevals + 1

      if (y1>=y.slice)
        break

      # Otherwise, shrink the approximation towards the current point.

      L[wt<0] <- wt[wt<0]
      U[wt>0] <- wt[wt>0]
    }

    # Save the accepted proposal.

    X[obs,] <- x1
    x0 <- x1
  }
  return(list(X=X, evals=nevals, S.eig=S.eig))
}

attr(oblique.hyperrect.sample, 'name') <- 'Oblique Hyperrect'

cheat.oblique.hyperrect.sample <- function(target.dist, x0, sample.size,
                                           tuning=1) {
  oblique.hyperrect.sample(target.dist, x0, sample.size, tuning, cheat=TRUE)
}

attr(cheat.oblique.hyperrect.sample, 'name') <- 'Cheat Oblique Hyperrect'
