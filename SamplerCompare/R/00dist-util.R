# From SamplerCompare, (c) 2010 Madeleine Thompson

# This file contains code for managing objects of the "scdist" class,
# which represent probability distributions in the SamplerCompare
# package.  The R help for these functions and the vignette "R/C Glue
# in SamplerCompare (doc/glue.pdf) may also be of interest.

# See ?make.dist for more information.

make.dist <- function(ndim, name, name.expression=NULL,
                      log.density=NULL, grad.log.density=NULL,
                      log.density.and.grad=NULL, initial=NULL,
                      mean=NULL, cov=NULL, mean.log.dens=NULL) {
  stopifnot(ndim>0)
  if (!is.null(mean))
    stopifnot(length(mean)==ndim)
  if (!is.null(cov))
    stopifnot(all(dim(cov)==ndim))
  stopifnot(is.null(log.density) || is.function(log.density))
  stopifnot(is.null(grad.log.density) || is.function(grad.log.density))
  stopifnot(is.null(log.density.and.grad) || is.function(log.density.and.grad))
  stopifnot(is.null(initial) || is.function(initial))

  # Define a scdist object as a list, filling in whatever the user gave us.

  ds <- list(ndim=ndim, name=name, name.expression=name.expression,
             log.density=log.density, grad.log.density=grad.log.density,
             log.density.and.grad=log.density.and.grad, initial=initial,
             mean=mean, cov=cov, mean.log.dens=mean.log.dens)

  # If the user gave us log.density.and.grad but not log.density,
  # fake up a log.density.

  if (is.null(log.density) && !is.null(log.density.and.grad)) {
    ds$log.density <- function(x) log.density.and.grad(x,FALSE)$log.density
  }

  # If the user gave us log.density.and.grad but not grad.log.density,
  # fake up a grad.log.density.

  if (is.null(grad.log.density) && !is.null(log.density.and.grad)) {
    ds$grad.log.density <-
      function(x) log.density.and.grad(x,TRUE)$grad.log.density
  }

  # If the user gave us log.density but not log.density.and.grad,
  # fake up a log.density.and.grad.  This will fail if compute.grad
  # is TRUE and grad.log.density was unspecified.

  if (is.null(log.density.and.grad) && !is.null(log.density)) {
    ds$log.density.and.grad <- function(x, compute.grad=FALSE) {
      if (compute.grad) {
        stopifnot(!is.null(ds$grad.log.density))
        return(list(log.density=ds$log.density(x),
                    grad.log.density=ds$grad.log.density(x)))
      } else {
        return(list(log.density=ds$log.density(x)))
      }
    }
  }

  # Mark the distribution with its class and return it.

  class(ds) <- 'scdist'
  return(ds)
}

# Overrides print() for scdist objects.  Prints basic information
# about a distribution.

print.scdist <- function(x, ...) {
  if (!is.null(x$c.log.density.and.grad))
    feat <- 'log-density and gradient implemented in C'
  else if (!is.null(x$log.density) && !is.null(x$grad.log.density))
    feat <- 'known log-density and gradient'
  else if (!is.null(x$log.density))
    feat <- 'known log-density'
  else
    feat <- 'unknown log-density'
  cat(sprintf("%s (%d dimensions, %s)\n", x$name, x$ndim, feat))
}

# See ?make.c.dist for more information.

make.c.dist <- function(ndim, name, c.log.density, c.context=NULL,
                        name.expression=NULL, mean=NULL, cov=NULL) {

  # Make a distribution object.

  ds <- make.dist(ndim, name, name.expression=name.expression,
                  mean=mean, cov=cov)
  
  # Locate the symbol for the log density.

  ndim <- as.integer(ndim)
  stopifnot(is.character(c.log.density) && length(c.log.density==1))
  ds$sym <- raw.symbol(c.log.density)

  # Define log.density.and.grad wrapper that calls the C version.

  ds$log.density.and.grad <- function(x, compute.grad=FALSE) {
    stopifnot(is.numeric(x) && length(x)==ndim)
    r <- .Call(R_invoked_C_glue, ds$sym, c.context, x, compute.grad)
    return(r)
  }

  # Define log.density and grad.log.density as wrappers.

  ds$log.density <- function(x) {
    ds$log.density.and.grad(x)$log.density
  }

  ds$grad.log.density <- function(x) {
    ds$log.density.and.grad(x, compute.grad=TRUE)$grad.log.density
  }

  # Note c.log.density.and.grad and c.context for debugging and print().

  ds$c.log.density.and.grad <- c.log.density
  ds$c.context <- c.context

  return(ds)
}

# Ensures that the log density and its gradient in a distribution
# object are consistent with each other at point x.  Uses the central
# difference formula on each coordinate in turn.
#
# See ?check.dist.gradient for more information.  ?make.dist also has
# an example.  See Nocedal & Wright, Numerical Optimization 2nd ed.
# sec. 8.1 (p. 196-197) for more information on the central difference
# formula.

check.dist.gradient <- function(ds, x, h=1e-7) {

  # Compute gradient with analytic formula.

  g <- ds$grad.log.density(x)
  stopifnot(all(is.finite(g)))
  stopifnot(length(g)==length(x))

  # Check each coordinate.

  for (i in 1:length(x))
    check.discrete.gradient(ds$log.density, x, g, i, h)
}

# Check analytic gradient of f, x.df, against numerical derivative
# in dimension i at point x.

check.discrete.gradient <- function(f, x, g, i, h) {

  # element vector

  e <- c(rep(0,i-1),1,rep(0,length(x)-i))

  # approximation

  g.discrete <- (f(x+e*h*x[i]) - f(x-e*h*x[i])) / 2 / h / x[i]
  stopifnot(is.finite(g.discrete))

  # If relative error is more than 0.001, fail.

  if((g[i]==0 && abs(g.discrete) > 1e-8) ||
     (g[i]!=0 && abs( (g.discrete-g[i])/g[i] ) > 1e-3)) {
    stop("Gradients do not match.\n",
         sprintf("analytic grad[%d] = %.5g   discrete grad[%d] = %.5g\n",
                 i, g[i], i, g.discrete))
  }
}

# Turn a list of step functions into a sampler function.

compounded.sampler <- function(step.functions, name, name.expr=NULL) {
  sampler <- function(target.dist, x0, sample.size, limit=NULL, ...) {
    # Initialization

    X <- array(NA,c(sample.size,length(x0)))
    step <- list(x=x0, y=target.dist$log.density(x0))
    grads <- 0
    evals <- 1

    # Loop that generates sample.size observations from target distribution

    for (obs in 1:sample.size) {
      
      # Take one step with each method.

      for (s in 1:length(step.functions)) {
        step <- step.functions[[s]](target.dist, step$x, step$y, ...)
        evals <- evals + step$evals
        grads <- grads + step$grads
      }

      # Store accepted observation and check to make sure the procedure
      # has not been running too long.

      X[obs,] <- step$x

      if (!is.null(limit) && evals > limit * sample.size) {
        X <- X[1:obs,]
        break
      }
    }
    return(list(X=X, evals=evals, grads=grads))
  }

  # Fill in attributes of sampler and return it.

  attr(sampler, 'name') <- name
  if (!is.null(name.expr))
    attr(sampler, 'name.expr') <- name.expr
  return(sampler)
}
