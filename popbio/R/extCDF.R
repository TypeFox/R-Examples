extCDF <- function(mu,sig2,Nc, Ne,tmax=50)
{
   ## some functions from box 3.3 (erf is matlab error function)
   erf <- function(y) 2 * pnorm(y * sqrt(2)) - 1
   phi <- function(z) 0.5 * (1 + (erf(z / sqrt(2))))
   # log distance from the quasi-extinction threshold
   d <- log(Nc/Ne) 
   G<-numeric(tmax)
   for (x in 1:tmax)
   {
      G[x] <- phi((-d-mu*x)/sqrt(sig2*x)) +
      exp(-2*mu*d/sig2)* phi((-d+mu*x)/sqrt(sig2*x))
   }
   G
}

 
