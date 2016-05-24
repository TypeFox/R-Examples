# 23.12.2009 version 1, mvdl
qqWeibullLimit <- function(y, p , iLambda, alpha)
{

   par  <- fitWeibull(y[iLambda],p[iLambda])
   yHat <- qweibull(p, par$k, par$lambda)
   res  <- log(y) - log(yHat)
   fac <- length(iLambda)/(length(iLambda)-2)
   sigmaE <- sqrt(fac*mean(res[iLambda]^2))
   
   L <- getLplusLmin(sigmaE, alpha)

   return(list(limit=c(Left=L$Lmin,Right=L$Lplus),
               residuals=res,
               sigmaE=sigmaE,
               k=par$k,
               lambda=par$lambda,
               R2=par$R2))
}
