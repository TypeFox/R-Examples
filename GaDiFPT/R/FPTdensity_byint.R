FPTdensity_byint <- 
  function(obj,n1max) 
    # obj is an inputlist class object yielding all the input parameters    
{
# General function for numerical integration of the Volterra integral equation
# with user-provided functions for infinitesimal drift and variance
  
# evaluation of FPT pdf and FPT cdf on timepoints up to Nfin 

   ntime <- floor((obj$Tfin - obj$t0)/obj$deltat)  
   if (n1max < ntime)  ntime <- n1max
   time.points <- seq(obj$t0,obj$Tfin,obj$deltat)
    
    # FPT density vector
   g0 <- numeric(ntime+1)
    
    # FPT distribution vector
   gg0 <- numeric(ntime+1)

   g0[1] <- 0
   for (i in 2:(ntime+1)) {
     g0[i] <- - psi(time.points[i],obj$x0,obj$t0);
     ti <- time.points[1:(i-1)]  
     gcurr <- psi(time.points[i],S(ti),ti)
     g0[i] <- g0[i] + sum(gcurr*g0[1:(i-1)])*obj$deltat
   }

  dt2 <- 0.5*obj$deltat

  gg0[1] <-  g0[1]*dt2
  gg0[2:(ntime+1)]  <- (g0[1:ntime]+g0[2:(ntime+1)])*dt2  
  gg0 <- cumsum(gg0)

  answer <- list(time.points,g0,gg0)
  attr(answer,"names") <- c("time","g0","gg0") 
  class(answer) <- c("FPTdensity", "list")
  return(answer)
}
