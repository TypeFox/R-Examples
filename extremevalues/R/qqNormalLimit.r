# 23.12.2009 version 1, mvdl
# 14.05.2010 solved bug in parameter passing, mvdl
qqNormalLimit <- function(y, p , iLambda, alpha)
{

   par  <- fitNormal(y[iLambda], p[iLambda])
   yHat <- qnorm(p, par$mu, par$sigma)
   res  <- y - yHat
   fac <- length(iLambda)/(length(iLambda)-2)
   sigmaE <- sqrt(fac*mean(res[iLambda]^2))
   
   L <- getLplusLmin(sigmaE, alpha)

   return(list(limit=c(Left=L$Lmin,Right=L$Lplus),
               residuals=res,
               sigmaE=sigmaE,
               mu=par$mu,
               sigma=par$sigma,
               R2=par$R2))
}
