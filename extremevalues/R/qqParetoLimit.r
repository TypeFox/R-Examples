# 23.12.2009 version 1, mvdl
qqParetoLimit <- function(y, p , iLambda, alpha)
{

   par  <- fitPareto(y[iLambda], p[iLambda])
   yHat <- qpareto(p, par$ym, par$alpha)
   res  <- log(y) - log(yHat)
   fac <- length(iLambda)/(length(iLambda)-2)
   sigmaE <- sqrt(fac*mean(res[iLambda]^2))
   
   L <- getLplusLmin(sigmaE, alpha)

   return(list(limit=c(Left=L$Lmin,Right=L$Lplus),
               residuals=res,
               sigmaE=sigmaE,
               ym=par$ym,
               alpha=par$alpha,
               R2=par$R2))
}
