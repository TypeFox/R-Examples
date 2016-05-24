# 19.10.2009, changed p-estimator
# 24.11.2009, added Weibull distribution
# 22.12.2009, added left limit, Changed rho default, added input checks.
# 08.01.2010, switched to Makkonen's equation for plot positions.
getOutliersI <- function(y, rho=c(1,1), FLim=c(0.1,0.9), distribution="normal")
{

   if ( !is.vector(y) ) 
      stop("First argument is not of type vector")
   if ( sum(y < 0) > 0 & !(distribution == "normal") )
      stop("First argument contains nonpositive values")
   if ( sum( rho <= 0, na.rm=TRUE ) > 0 )
      stop("Values of rho must be positive")
   if ( FLim[2] <= FLim[1] | sum( FLim < 0 | FLim > 1) >0 )
      stop("Invalid range in FLim: 0<=FLim[1]<FLim[2]<=1")
   if ( ! distribution %in% c("lognormal", "pareto", "exponential", "weibull", "normal") )
      stop("Invalid distribution (lognormal, pareto, exponential, weibull, normal).")

   Y <- y;
 
   y <- sort(y);
   N <- length(y)
   P <- (1:N)/(N+1)
   Lambda <- P >= FLim[1] & P<=FLim[2]
 
   y <- y[Lambda];
   p <- P[Lambda];
   out <- switch(distribution,
         lognormal = getLognormalLimit(y, p, N, rho),
         pareto = getParetoLimit(y, p, N, rho),
         exponential = getExponentialLimit(y, p, N, rho),
         weibull = getWeibullLimit(y, p, N, rho),
         normal = getNormalLimit(y, p, N, rho)
         )
   
   out$method <- "Method I"
   out$distribution=distribution
   out$iRight = which( Y > out$limit[2] )
   out$iLeft = which( Y < out$limit[1] )
   out$nOut = c(Left=length(out$iLeft), Right=length(out$iRight))
   out$yMin <- y[1]
   out$yMax <- tail(y,1)
   out$rho = c(Left=rho[1], Right=rho[2])
   out$Fmin = FLim[1]
   out$Fmax = FLim[2]

   return(out);
}

