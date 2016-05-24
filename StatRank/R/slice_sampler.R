################################################################################
# Slice Sampling
################################################################################


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

# uni.slice.calls <- 0	# Number of calls of the slice sampling function
# uni.slice.evals <- 0	# Number of density evaluations done in these calls


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

uni.slice <- function (x0, g, w=1, m=Inf, lower=-Inf, upper=+Inf, gx0=NULL)
{
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
  { #uni.slice.evals <<- uni.slice.evals + 1
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


# FUNCTION TO TEST THE UNI.SLICE FUNCTION.  
#
# Produces Postscript plots in slice-test.ps, one page per test, with the
# tests described in the code below.  
#
# Each test applies a series of univariate slice sampling updates (ss*thin of 
# them) to some distribution, starting at a point drawn from that distribution,
# with particular settings of the slice sampling options.  The page for a 
# test contains the following:
# 
#   - a trace plot of the results (at every 'thin' updates)
#   - a plot of the autocorrelations for this trace
#   - a plot of the bivariate distribution before and after 'thin' updates
#   - a qqplot of the sample produced vs. a correct sample
#   - the average number of evaluations per call
#   - the result of a t test for the sample mean vs. the correct mean, based
#     on 200 equally spaced points from the sample generated (which are
#     presumed to be virtually independent)

# uni.slice.test <- function ()
# {
#   postscript("slice-test.ps")
#   par(mfrow=c(2,2))

#   # Function to do the slice sampling updates.  

#   updates <- function (x0, g, reuse=FALSE)
#   { 
#     uni.slice.calls <<- 0
#     uni.slice.evals <<- 0

#     s <<- numeric(ss)
#     x1 <- x0
#     s[1] <<- x0
#     last.g <- NULL

#     for (i in 2:ss)
#     { for (j in 1:thin)
#       { if (reuse)
#         { x1 <- uni.slice (x1, g, w=w, m=m, lower=lower, upper=upper, 
#                            gx0=last.g)
#           last.g <- attr(x1,"log.density")
#         }
#         else
#         { x1 <- uni.slice (x1, g, w=w, m=m, lower=lower, upper=upper)
#         }
#       }
#       s[i] <<- x1
#     }
#   }

#   # Function to display the results.

#   display <- function (r,mu,test)
#   { 
#     plot(s,type="p",xlab="Iteration",ylab="State",pch=20)

#     title (paste( test, "  ss =",ss," thin =",thin))

#     acf (s, lag.max=length(s)/20, main="")

#     title (paste ("w =",w," m =",m," lower =",lower," upper =",upper))

#     plot(s[-1],s[-length(s)],pch=20,xlab="Current state",ylab="Next state")

#     title (paste ("Average number of evaluations:",
#                    round(uni.slice.evals/uni.slice.calls,2)))

#     qqplot(r,s,pch=".",
#        xlab="Quantiles from correct sample",
#        ylab="Quantiles from slice sampling")
#     abline(0,1)

#     p.value <- t.test (s[seq(1,length(s),length=200)]-mu) $ p.value
#     title (paste ("P-value from t test:",round(p.value,3)))
#   }

#   # Standard normal, m = Inf.

#   set.seed(1)

#   ss <- 2000
#   thin <- 3
#   w <- 1.5
#   m <- Inf
#   lower <- -Inf
#   upper <- +Inf

#   updates (rnorm(1), function (x) -x^2/2)
#   display (rnorm(ss),0,"Standard normal")

#   # Standard normal, reusing density, m = Inf.

#   set.seed(1)

#   ss <- 2000
#   thin <- 3
#   w <- 1.5
#   m <- Inf
#   lower <- -Inf
#   upper <- +Inf

#   updates (rnorm(1), function (x) -x^2/2, reuse=TRUE)
#   display (rnorm(ss),0,"Standard normal, reusing density")

#   # Normal mixture, m = 1.

#   set.seed(1)

#   ss <- 2000
#   thin <- 3
#   w <- 2.2
#   m <- 1
#   lower <- -Inf
#   upper <- +Inf

#   updates (rnorm(1,-1,1), function (x) log(dnorm(x,-1,1)+dnorm(x,1,0.5)))
#   display (c (rnorm(floor(ss/2),-1,1), rnorm(ceiling(ss/2),1,0.5)), 0,
#            "Normal mixture")

#   # Normal mixture, m = 3.

#   set.seed(1)

#   ss <- 2000
#   thin <- 3
#   w <- 1.8
#   m <- 3
#   lower <- -Inf
#   upper <- +Inf

#   updates (rnorm(1,-1,1), function (x) log(dnorm(x,-1,1)+dnorm(x,1,0.5)))
#   display (c (rnorm(floor(ss/2),-1,1), rnorm(ceiling(ss/2),1,0.5)), 0,
#            "Normal mixture")

#   # Exponential, m = Inf.

#   set.seed(1)

#   ss <- 2000
#   thin <- 3
#   w <- 10
#   m <- Inf
#   lower <- 0
#   upper <- +Inf

#   updates (rexp(1), function (x) -x)
#   display (rexp(ss), 1, "Exponential")

#   # Exponential, m = 2.

#   set.seed(1)

#   ss <- 2000
#   thin <- 3
#   w <- 1.5
#   m <- 2
#   lower <- 0
#   upper <- +Inf

#   updates (rexp(1), function (x) -x)
#   display (rexp(ss), 1, "Exponential")

#   # Beta(0.5,0.8).

#   set.seed(1)

#   ss <- 2000
#   thin <- 3
#   w <- 1e9
#   m <- Inf
#   lower <- 0
#   upper <- 1

#   updates (rbeta(1,0.5,0.8), function (x) dbeta(x,0.5,0.8,log=TRUE))
#   display (rbeta(ss,0.5,0.8), 0.5/(0.5+0.8), "Beta(0.5,0.8)")

#   dev.off()
# }
