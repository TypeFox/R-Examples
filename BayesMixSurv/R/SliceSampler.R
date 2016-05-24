# wrapper around a multivariate function to convert it to univariate, to be used with slice sampler
bayesmixsurv.convert.multivar.to.univar <- function(xk, k, x, f, ...) {
  x[k] <- xk
  return (f(x, ...))
}
# Gibbs wrapper around univariate slice sampler to convert it to a multivariate sampler
bayesmixsurv.multislice.from.unislice <- function(xm, fm, ..., w=1.0, m=Inf, lower=rep(-Inf,length(xm)), upper=rep(+Inf,length(xm))) {
  K <- length(xm)
  for (k in 1:K) {
    xm[k] <- bayesmixsurv.uni.slice(xm[k], bayesmixsurv.convert.multivar.to.univar, k, xm, fm, ..., w=w, m=m, lower=lower[k], upper=upper[k])
  }
  return (xm)
}

# R FUNCTIONS FOR PERFORMING UNIVARIATE SLICE SAMPLING.
#
# Radford M. Neal, 17 March 2008.
#
# Implements, with slight modifications and extensions, the algorithm described
# in Figures 3 and 5 of the following paper:
#
#   Neal, R. M (2003) "Slice sampling" (with discussion), Annals of Statistics,
#      vol. 31, no. 3, pp. 705-767.
#
# See the documentation for the function uni.slice below for how to use it.
# The function uni.slice.test was used to test the uni.slice function.


# GLOBAL VARIABLES FOR RECORDING PERFORMANCE.

#uni.slice.calls <- 0  # Number of calls of the slice sampling function
#uni.slice.evals <- 0	# Number of density evaluations done in these calls


# UNIVARIATE SLICE SAMPLING WITH STEPPING OUT AND SHRINKAGE.
#
# Performs a slice sampling update from an initial point to a new point that 
# leaves invariant the distribution with the specified log density function.
#
# Arguments:
#
#   x0    Initial point
#   g     Function returning the log of the probability density (plus constant)
#   w     Size of the steps for creating interval (default 1)
#   m     Limit on steps (default infinite)
#   lower Lower bound on support of the distribution (default -Inf)
#   upper Upper bound on support of the distribution (default +Inf)
#   gx0   Value of g(x0), if known (default is not known)
#
# The log density function may return -Inf for points outside the support 
# of the distribution.  If a lower and/or upper bound is specified for the
# support, the log density function will not be called outside such limits.
#
# The value of this function is the new point sampled, with an attribute
# of "log.density" giving the value of the log density function, g, at this
# point.  Depending on the context, this log density might be passed as the 
# gx0 argument of a future call of uni.slice. 
#
# The global variable uni.slice.calls is incremented by one for each call
# of uni.slice.  The global variable uni.slice.evals is incremented by the
# number of calls made to the g function passed.
#
# WARNING:  If you provide a value for g(x0), it must of course be correct!
# In addition to giving wrong answers, wrong values for gx0 may result in
# the uni.slice function going into an infinite loop.

bayesmixsurv.uni.slice <- function (x0, f, ..., w=1, m=0, lower=-Inf, upper=+Inf, gx0=NULL)
{
  g <- function(x) f(x,...)
  # Check the validity of the arguments.
  
  if (!is.numeric(x0) || length(x0)!=1
      || !is.function(g) 
      || !is.numeric(w) || length(w)!=1 || w<=0 
      || !is.numeric(m) || !is.infinite(m) && (m<=0 || m>1e9 || floor(m)!=m)
      || !is.numeric(lower) || length(lower)!=1 || x0<lower
      || !is.numeric(upper) || length(upper)!=1 || x0>upper
      || upper<=lower 
      || !is.null(gx0) && (!is.numeric(gx0) || length(gx0)!=1))
  { 
    stop ("Invalid slice sampling argument")
  }
  
  # Keep track of the number of calls made to this function.
  
  #uni.slice.calls <<- uni.slice.calls + 1
  
  # Find the log density at the initial point, if not already known.
  
  if (is.null(gx0)) 
  {
    #uni.slice.evals <<- uni.slice.evals + 1
    gx0 <- g(x0)
  }
  
  # Determine the slice level, in log terms.
  
  logy <- gx0 - rexp(1)
  
  # Find the initial interval to sample from.
  
  u <- runif(1,0,w)
  L <- x0 - u
  R <- x0 + (w-u)  # should guarantee that x0 is in [L,R], even with roundoff
  
  # Expand the interval until its ends are outside the slice, or until
  # the limit on steps is reached.
  
  if (is.infinite(m))  # no limit on number of steps
  { 
    repeat
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
    }
    
    repeat
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
    }
  }
  
  else if (m>1)  # limit on steps, bigger than one
  { 
    J <- floor(runif(1,0,m))
    K <- (m-1) - J
    
    while (J>0)
    { if (L<=lower) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(L)<=logy) break
      L <- L - w
      J <- J - 1
    }
    
    while (K>0)
    { if (R>=upper) break
      #uni.slice.evals <<- uni.slice.evals + 1
      if (g(R)<=logy) break
      R <- R + w
      K <- K - 1
    }
  }
  
  # Shrink interval to lower and upper bounds.
  
  if (L<lower) 
  { L <- lower
  }
  if (R>upper)
  { R <- upper
  }
  
  # Sample from the interval, shrinking it on each rejection.
  
  repeat
  { 
    x1 <- runif(1,L,R)
    
    #uni.slice.evals <<- uni.slice.evals + 1
    gx1 <- g(x1)
    
    if (gx1>=logy) break
    
    if (x1>x0) 
    { R <- x1
    }
    else 
    { L <- x1
    }
  }
  
  # Return the point sampled, with its log density attached as an attribute.
  
  attr(x1,"log.density") <- gx1
  return (x1)
  
}
