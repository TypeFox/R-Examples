# file  : getExpLimit.r
# author: Mark van der Loo (mark.vanderloo@gmail.com)
# 
# Determine outlier limit based exponential distribution.
#
# INPUT
# y     : vector of observed values between pmin and pmax
# p     : vector of observed quantiles (y_i estimates the p_i'th quantile)
# N     : total number of observations
# rho   : outlier parameters
#
# OUTPUT (list)
# lambda: estimate parameter
# R2    : R-squared value of fit.
#
# History
# 22.10.2009    version 1
# 22.12.2009    version 2 (mvdl) added left limit.
# 


getExponentialLimit <- function(y, p, N, rho)
{
   param <- fitExponential(y, p)
   ell <- c(Left=0, Right=Inf)
   if ( !is.na(rho[1]) )
    ell[1] <- -log(1-rho[1]/N)/param$lambda
   if ( !is.na(rho[2]) )
    ell[2] <- log(N/rho[2])/param$lambda

   return(list(lambda=param$lambda, 
               R2=param$R2,
               nFit=length(y),
               limit=ell)
         )
}
