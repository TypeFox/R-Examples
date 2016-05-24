FPTsimul <- 
  function(obj,M) 
    {
# Implementation of the Hazard Rate Method
# with user-provided functions for the drift, infinitesimal variance and threshold
    spikes <- numeric(M)
    
    time.points <- obj$time
    rate.vec <- obj$g0/(1-obj$gg0)
    N1time <- length(obj$time)
    t0 <- time.points[1]
    deltat <- time.points[2]-time.points[1]
    
# refines N1 to Nmax (start of stationary regime) to shorten computations
    
   Nmin <- floor(N1time/10)
   Nmax <- which.min(abs(rate.vec[(Nmin+1):(N1time+1)]-rate.vec[Nmin:N1time]))+Nmin
                  

# fixed tolerance

   tol <- 0.0001

   tmax <- t0 + Nmax*deltat
   lamax <- 1.01*max(rate.vec)
   nreject <- 0
   for (ifire in 1:M) {
     Tcurr <- t0
     repeat {
     Tcurr <- Tcurr - log(runif(1))/lamax
     lambdacurr <- rate.vec[Nmax]
     if (Tcurr < tmax) {
       icurr <- min(floor((Tcurr-t0)/deltat+1.001),(N1time-1))
       lambdacurr <- rate.vec[icurr]+(rate.vec[icurr+1]-rate.vec[icurr])*(Tcurr-time.points[icurr])/deltat
      }
      u2 <- runif(1)  
      lacheck <- abs(lambdacurr/lamax - u2)
      if (lacheck <= tol ) {
      nreject <- nreject + 1
      Tcurr <- t0 
      }
      if ((lacheck > tol) & (lambdacurr/lamax >= u2) ) break()
  }
  
  spikes[ifire] <- Tcurr
}

return(spikes)
}
