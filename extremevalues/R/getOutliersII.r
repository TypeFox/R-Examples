# 24.12.2009 version 1, mvdl
# 08.01.2010 changed to Makkonen's solution for plot positions
getOutliersII <- function(y, alpha=c(0.05, 0.05), FLim=c(0.1, 0.9), distribution="normal", returnResiduals=TRUE)
{
# Input check
   if ( !is.vector(y) ) 
      stop("First argument is not of type vector")
   if ( sum(y < 0) > 0 & !(distribution == "normal") )
      stop("First argument contains nonpositive values")
   if ( sum( alpha <= 0 | alpha >= 1, na.rm=TRUE ) > 0 )
      stop("Values of alpha must be between 0 and 1")
   if ( FLim[2] <= FLim[1] | sum( FLim < 0 | FLim > 1) >0 )
      stop("Invalid range in FLim: 0<=FLim[1]<FLim[2]<=1")
   if ( ! distribution %in% c("lognormal", "pareto", "exponential", "weibull", "normal") )
      stop("Invalid distribution (lognormal, pareto, exponential, weibull, normal).")

# prepare for regression
   iy <- order(y) 
   Y <- y[iy]
   N <- length(Y)
   p <- seq(1,N)/(N+1)
   iLambda <- which( p >= FLim[1] & p <= FLim[2] )
   if (length(iLambda) <= 2 )
      stop("Number of observations in fit is too small to estimate variance of residuals (need at least 3)")
# regression and limit calculation
   par <- switch(distribution,
      lognormal   = qqLognormalLimit(Y, p, iLambda, alpha),
      weibull     = qqWeibullLimit(Y, p, iLambda, alpha),
      pareto      = qqParetoLimit(Y, p, iLambda, alpha),
      exponential = qqExponentialLimit(Y, p, iLambda, alpha),
      normal      = qqNormalLimit(Y, p, iLambda, alpha),
      )

# locate outliers
   iLmin <- iLplus <- numeric(0)
   i <- N + 1
   while ( par$residuals[i-1] > par$limit[2] & i > tail(iLambda,1)+1 )
      i <- i-1
   if ( i <= N )
      iLplus <- iy[i:N]

   i <- 0
   while ( par$residuals[i+1] < par$limit[1] & i < iLambda[1]-1 )
      i <- i+1
   if ( i > 0 )
      iLmin <- iy[1:i]

# organize output
   out <- par
   out$method <- "Method II"
   out$distribution <- distribution
   out$iRight <- iLplus
   out$iLeft <- iLmin 
   out$nOut <- c(Left=length(iLmin),Right=length(iLplus))
   out$yMin <- head(Y[iLambda],1)
   out$yMax <- tail(Y[iLambda],1)
   out$alphaConf<-c(Left=alpha[1], Right=alpha[2])
   out$nFit <- length(iLambda)
   if ( returnResiduals ){
      out$residuals <- numeric(N)
      out$residuals[iy] <- par$residuals
      }
   return(out)
}





