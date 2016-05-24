# 23.12.2009 version 1, mvdl
qqExponentialLimit <- function(y, p , iLambda, alpha)
{

   par  <- fitExponential(y[iLambda], p[iLambda])
   yHat <- (-1/par$lambda)*log(1-p)
   res  <- y - yHat
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
