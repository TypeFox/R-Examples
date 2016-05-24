# file  : getWeibulLimit.r
# author: Mark van der Loo (mark.vanderloo@gmail.com)
# 
# Determine outlier limit assuming Weibull distribution.
#
# INPUT
# y     : vector of observed values between pmin and pmax
# p     : vector of observed quantiles (y_i estimates the p_i'th quantile)
# N     : total number of observations
# rho   : outlier parameter
#
# OUTPUT (list)
# limit : estimate parameter
# R2    : R-squared value of fit.
#
# History
# 24.11.2009    version 1
# 22.12.2009    version 2 (mvdl) added left limit.
#

getWeibullLimit <- function(y, p, N, rho)
{
   param <- fitWeibull(y,p)
   ell <- c(Left=0, Right=Inf)
   if ( !is.na(rho[1]) )
      ell[1] <- param$lambda*(-log(1-rho[1]/N))^(1/param$k)
   if ( !is.na(rho[2]) )
      ell[2] <- param$lambda*log(N/rho[2])^(1/param$k)
   return(list(k=param$k, 
               lambda=param$lambda,
               nFit=length(y),
               R2=param$R2,
               limit=ell)
         )
}
