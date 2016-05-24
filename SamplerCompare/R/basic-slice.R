# From SamplerCompare, (c) 2010 Madeleine Thompson

# This file contains implementations of several slice samplers from:
#
#   Radford Neal (2003), Slice Sampling, Annals of Statistics 31(3):705-767.

# hyperrectangle.sample implements the method of Neal (2003) sec. 5.1.
# The R help page has some more information.

hyperrectangle.sample <- function(target.dist, x0, sample.size, tuning=1,
                             use.gradient=TRUE, limit=length(x0)*100) {
  X <- array(NA,c(sample.size,target.dist$ndim))
  X[1,] <- x0
  nevals <- 0
  ngrads <- 0

  # Use undocumented API for distributions to mark their lower
  # bounds explicitly.

  if (!is.null(target.dist$lower.bound)) {
    lower.bound <- target.dist$lower.bound
  } else {
    lower.bound <- rep(-Inf, target.dist$ndim)
  }

  for (obs in 2:sample.size) {

    # Set x0 to the last state and y.slice to a new slice level.
    # Make lower and upper the bounds of a randomly positioned box
    # around x0.

    x0 <- X[obs-1,,drop=TRUE]
    nevals <- nevals + 1
    y.slice <- target.dist$log.density(x0) - rexp(1)
    lower <- x0 - runif(target.dist$ndim) * tuning
    upper <- lower + tuning
    stopifnot(all(upper>lower.bound))
    lower <- pmax(lower, lower.bound)

    repeat {

      # Draw a proposal, x1, from the box bounded by lower and
      # upper.  If the proposal is inside the slice, accept it.

      x1 <- lower + (upper-lower) * runif(target.dist$ndim)
      nevals <- nevals + 1
      if (target.dist$log.density(x1)>=y.slice)
        break

      # If we are using gradients, shrink the box in the direction
      # that the gradient implies the biggest change in log density
      # across that direction.  If we are not using gradients, shrink
      # in all directions.

      if (use.gradient) {
        ngrads <- ngrads + 1
        grad <- array(target.dist$grad.log.density(x1), c(target.dist$ndim,1))
        # FIXME: This variant needs justification, e.g. funnel with high tuning
        shrink.indices <- which.max(abs(grad)*(upper-lower))
      } else {
        shrink.indices <- 1:target.dist$ndim
      }

      # Shrink the box toward x0 in each direction we are changing it.

      for (i in shrink.indices) {
        if (x1[i] > x0[i])
          upper[i] <- x1[i]
        else
          lower[i] <- x1[i]
      }

      if (nevals > limit * sample.size)
        break
    }

    # Store the accepted proposal and make sure we have not run too long.

    X[obs,] <- x1

    if (nevals > limit * sample.size) {
      X <- X[1:(obs-1),,drop=FALSE]
      break
    }
  }
  return(list(X=X, evals=nevals, grads=ngrads))
}

attr(hyperrectangle.sample, 'name') <- 'Hyperrectangle'


# nograd.hyperrectangle.sample is a wrapper around hyperrectangle.sample
# with use.gradient set to false so that it can be passed directly
# to compare.samplers.

nograd.hyperrectangle.sample <- function(...) {
  return(hyperrectangle.sample(..., use.gradient=FALSE))
}

attr(nograd.hyperrectangle.sample, 'name') <- 'Hyperrectangle (no grad)'


# stepout.slice.sample implements the stepping-out procedure of
# Neal (2003) sec. 4.1.  The R help page has some more information.
# If step.out is set to FALSE, no stepping out is performed and an
# interval of fixed size is used as an initial estimate of the slice.

stepout.slice.sample <- function(target.dist, x0, sample.size, tuning=1,
                                 step.out=TRUE, limit=length(x0)*100) {
  p <- length(x0)
  X <- array(NA,c(sample.size,p))
  X[1,] <- x0
  nevals <- 0

  for (obs in 2:sample.size) {

    # Set x0 to current state and y.slice to a new slice level.

    x0 <- X[obs-1,,drop=TRUE]
    y.slice <- target.dist$log.density(x0) - rexp(1)
    nevals <- nevals + 1

    # Transition each coordinate in turn.

    for (i in 1:p) {

      # If step.out is set, find the boundaries of the slice in
      # this coordinate using stepout.one.coord.  Otherwise, choose
      # a random interval.

      if (step.out) {
	lr <- stepout.one.coord(target.dist$log.density, x0, i, y.slice,
	                        tuning, limit*sample.size - nevals)
        nevals <- nevals + lr$nevals
        if (is.null(lr$left))
          return(list(X=X[1:(obs-1),,drop=FALSE], evals=nevals))
      } else {
        lr <- list(left = x0[i] - runif(1) * tuning)
        lr$right <- lr$left + tuning
      }

      # Draw proposals, shrinking the slice estimate along the way,
      # until a proposal is accepted.

      repeat {

	# Draw a proposal by modifying the current coordinate of
	# x0 to be uniformly chosen from the current slice estimate..

        x1 <- x0
        x1[i] <- lr$left + (lr$right-lr$left) * runif(1)

	# If the proposal is in the slice, accept it and move on
	# to the next coordinate.

        y1 <- target.dist$log.density(x1)
        nevals <- nevals + 1
        if (y1>=y.slice) {
          x0 <- x1
          break
        }

        # Otherwise, shrink the slice towards x0.

        if (x1[i] < x0[i])
          lr$left <- x1[i]
        else
          lr$right <- x1[i]
      }
    }

    # Having updated each coordinate, save the accepted proposal.
  
    X[obs,] <- x1

    # Make sure we haven't run too long.

    if (nevals > limit * sample.size) {
      X <- X[1:obs,,drop=FALSE]
      break
    }
  }
  return(list(X=X, evals=nevals))
}

attr(stepout.slice.sample, 'name') <- 'Step-out Slice'

# interval.slice.sample is a wrapper around stepout.slice.sample
# with step.out set to FALSE so that it can be passed directly
# to compare.samplers.

interval.slice.sample <- function(...) {
  return(stepout.slice.sample(..., step.out=FALSE))
}

attr(interval.slice.sample, 'name') <- 'Interval Slice'


# stepout.one.coord performs the stepping-out procedure on a single
# coordinate, returning an estimate of the slice endpoints.  The
# arguments are:
#
#   L           log.density function
#   x           current point
#   i           index of coordinate to estimate interval in
#   y.slice     slice level
#   w           initial slice size estimate
#   max.evals   maximum number of times to call L before aborting

stepout.one.coord <- function(L, x, i, y.slice, w, max.evals=NULL) {

  # Choose left and right to be the same as x, but with one coordinate
  # changed in each of them so that together they bound a randomly
  # positioned interval around x in coordinate i.

  left <- x
  left[i] <- left[i] - runif(1) * w
  right <- x
  right[i] <- left[i] + w

  # Compute the log density at each end of the proposed interval.

  left.y <- L(left)
  right.y <- L(right)
  nevals <- 2

  # Keep expanding as long as either end has log density smaller
  # than the slice level.

  while (left.y > y.slice || right.y > y.slice) {

    # Expand each side of the proposed slice with equal probability.

    if (runif(1) > 0.5) {
      right <- right + w
      right.y <- L(right)
    } else {
      left <- left - w
      left.y <- L(left)
    }
    nevals <- nevals + 1

    # Make sure we haven't run too long.  If we have, return without
    # an interval.

    if (!is.null(max.evals) && nevals >= max.evals)
      return(list(nevals=nevals))
  }

  # We found an acceptable interval; return it.

  return(list(left=left[i], right=right[i], nevals=nevals))
}
