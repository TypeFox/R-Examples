# file  : getParetoLimit.r
# author: Mark van der Loo (mark.vanderloo@gmail.com)
# 
# Determine outlier limit assuming pareto distribution.
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
# 22.10.2009    version 1
# 22.12.2009    version 2 (mvdl) added left limit.
#


getParetoLimit <- function(y, p, N, rho)
{
   param <- fitPareto(y,p)
   ell <- c(Left=-Inf, Right=Inf)
   if ( !is.na(rho[1]) )
      ell[1] <- param$ym*(1-rho[2]/N)^{-1/param$alpha}
   if ( !is.na(rho[2]) )
      ell[2] <- param$ym*(N/rho[2])^{1/param$alpha}
   
   return(list(ym=param$ym, 
               alpha=param$alpha,
               nFit=length(y),
               R2=param$R2, 
               limit=ell)
         )
}
