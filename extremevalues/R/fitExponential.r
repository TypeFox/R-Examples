# file  : fitPareto.r
# author: Mark van der Loo (mark.vanderloo@gmail.com)
#
# Determine parameter lambda (scale and shape)
# of a exponential distributed variable Y, by fitting (part of) 
# the cdf to an observed cdf.
#
# INPUT
# y     : vector of observed values
# p     : vector of observed quantiles (y_i estimates the p_i'th quantile)
#
# OUTPUT (list)
# lambda: estimate parameter
# R2    : R-squared value of fit.
#
# History
# 22.10.2009    version 1
# 05.12.2009    fixed bug in fit of lambda

fitExponential <- function(y,p)
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


#   Lambda <- -sum(log(1-p))/sum(y);
   Lambda <- -sum(log(1-p)^2)/sum(y * log(1-p))
   r2 <- 1 - var(y - (-log(1-p)/Lambda) )/var(y);

   return(list(lambda=Lambda, R2=r2));


}



