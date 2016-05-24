# file  : fitPareto.r
# author: Mark van der Loo (mark.vanderloo@gmail.com)
#
# Determine parameters ym (scale) and alpha (shape)
# of a pareto distributed variable Y, by fitting (part of) 
# the log(cdf) to an observed log(cdf).
#
# INPUT
# y     : vector of observed values
# p     : vector of observed quantiles (y_i estimates the p_i'th quantile)
#
# OUTPUT (list)
# ym    : estimate of scale parameter
# alpha : estimate of shape parameter
# R2    : R-squared value of fit. (logarithmic)
#
# History
# 22.10.2009    version 1
#

fitPareto <- function(y,p)
{
   if ( !is.vector(y) ) 
      stop("First argument is not of type vector")
   if ( sum(y<=0) > 0 )
      stop("First argument contains nonpositive values")
   if ( !is.vector(p)) 
      stop("First argument is not of type vector")
   if ( sum(p<=0) > 0 | sum(p>=1) >0 )
      stop("Second argument contains values out of range (0,1)")
   if (length(y) != length(p))
      stop("First and second argument have different length");

   N <- length(y);
   lnY <- as.matrix(log(y),nrow=N)
   p <- as.matrix(p,nrow=N)

   A <- matrix(0,nrow=N,ncol=2)
   A[,1] <- 1 + double(N);
   A[,2] <- log(1-p);
   param <- solve(t(A) %*% A)%*%t(A)%*%lnY
   r2 <- 1 - var(exp(A%*%param) - y)/var(y);
   
   return(list(ym=exp(param[1]), alpha=-1/param[2], R2=r2));
}

