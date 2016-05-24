# file  : getLognormalLimit.r
# author: Mark van der Loo (mark.vanderloo@gmail.com)
# 
# Determine outlier limit assuming lognormal distribution.
#
# INPUT
# y     : vector of observed values between pmin and pmax
# p     : vector of observed quantiles (y_i estimates the p_i'th quantile)
# N     : total number of observations
# rho   : outlier parameter
#
# OUTPUT (list)
# lambda: estimate parameter
# R2    : R-squared value of fit.
#
# History
# 03.12.2009    version 1
# 22.12.2009    version 2 (mvdl) added left limit.
#

getNormalLimit <- function(y, p, N, rho)
{
   param <- fitNormal(y,p)
   ell <- c(Left=-Inf, Right=Inf)
   if ( !is.na(rho[1]) )
      ell[1] <- sqrt(2)*param$sigma*invErf(2*rho[1]/N-1)+param$mu
   if ( !is.na(rho[2]) )
      ell[2] <- sqrt(2)*param$sigma*invErf(1-2*rho[2]/N)+param$mu
   return(list(mu=param$mu, 
               sigma=param$sigma,
               nFit=length(y),
               R2=param$R2,
               limit=ell)
         )
}
