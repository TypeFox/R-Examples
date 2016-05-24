forecast <- function(skeleton,steps) {
  #*********************************************************************
  #  This function uses the GPR interpolator to interpolate the flow 
  #  field and provide an interative step-by-step forecast 
  #  Regression.
  #
  #  Input: skeleton - Matrix containing the data skeleton
  #         steps - Number of steps to forecast   
  #
  #  Output: fcast - forecast in incremental steps of size equal to the 
  #                  knot spacing of the PSR
  #
  #  References: 1. C. E. Rasmussen and C. K. I. Williams, Gaussian Processes
  #                 for Machine Learning, Cambridge, MA, MIT Press, 2006
  #
  #              2. Frey, MR and Caudle, KA “Flow field forecasting for 
  #                 univariate time series,” Statistical Analysis and Data 
  #                 Mining, 2013 
  #
  #**********************************************************************
  knots <- length(skeleton$sd)-4
  space <- skeleton[knots+2,2]
  lknot <- skeleton[knots+3,2]
  sigma2 <- skeleton[knots+4,1] # Variance estimate of PSR noise
  
  sd <- skeleton$sd[1:(knots+1)]
  delta <- skeleton$delta[1:(knots+1)]
  responses <- delta[4:(knots+1)]
  
  ssd = sd(sd)        # Standard deviation of the systematic components
  sdelta = sd(delta)  # Standard deviation of the slopes
  
  #**********************************************************************
  #  Create two matrices to hold the current and two previous slopes
  #  and one matrix to hold the previous three systematicallly determine
  #  components (SDCs)
  #
  #  hist.sd has size (knots - 2) x 3
  #  hist.delta has size (knots - 2 ) x 3
  #
  #  For example, the 10th row of hist.sd would be: 
  #                  
  #                   sd[11]  sd[12]  sd[13]
  #
  #               the 10th row of hist.delta would be:
  #
  #                   delta[10]  delta[11]  delta[12]
  #**********************************************************************
  
  # Holds the 3 previous SDCs
  hist.sd <-  matrix(data=0,nrow=knots-2,ncol=3)
  
  # Holds the current and 2 previous slopes
  hist.delta <- matrix(data=0,nrow=knots-2,ncol=3)
  
  # Holds the GPR response associated with the previous SDCs and slopes
  responses <- delta[4:(knots+1)]
  
  # Extract the history space
  for (i in 1:(knots-2)){
    hist.sd[i,1] <- skeleton$sd[i+1]
    hist.sd[i,2] <- skeleton$sd[i+2]
    hist.sd[i,3] <- skeleton$sd[i+3]
    hist.delta[i,1] <- skeleton$delta[i]
    hist.delta[i,2] <- skeleton$delta[i+1]
    hist.delta[i,3] <- skeleton$delta[i+2]
  }
  
  #******************** Prepare to Forecast ***************************
  
  tf <- rep(0,steps)        # Vector to hold forecast times
  fcast <- rep(0,steps)     # Vector to hold forecast values
  
  tf[1] <- lknot + space     
  fcast[1] = sd[knots] + space*delta[knots] 
  
  # Extract the current history
  rec3.sd <- c(sd[knots-1], sd[knots], fcast[1])         
  rec3.delta <- c(delta[knots-2], delta[knots-1], delta[knots]) 
  
  #******************** Flow Field Forecast ***************************
  
  std.error <- rep(0,steps)        # Vector to hold the standard errors 
  std.error[1] <-  sqrt(sigma2)
  
  h <- matrix(data=0,nrow=knots-2,ncol=6) 
  h <- cbind(hist.sd,hist.delta)
  
  omega <- rep(0,steps-1)
  k <-0
  for (i in (2:steps)){
    # Estimate GPR predicted change
    gpr.change <- gpr(h,rec3.sd,rec3.delta,ssd,sdelta,responses)
    kmat <- gpr.change$first
    kvec <- gpr.change$second
    ytilda <- solve(kmat,responses)
    new.delta = t(kvec) %*% ytilda
    responses <- c(responses,new.delta)
    omega[i-1] <- t(kvec)%*%solve(kmat,kvec)
    
    fcast[i] = fcast[i-1] + space*new.delta  # new predicted level
    row <- length(h[,1])
    
    # Newest 3 levels and slopes
    rec3.sd <- c(h[row,2], fcast[i-1], fcast[i])   
    rec3.delta <- c(h[row,5], h[row,6], new.delta) 
    newdata <- c(rec3.sd,rec3.delta)
    h <- rbind(h,newdata)
    
    tf[i] = tf[i-1]+ space    # updated future times
    k <- k + t(kvec)%*%solve(kmat,kvec)
    std.error[i] <- sqrt(sigma2 + (tf[i] - tf[1])*space*sdelta^2*(k/i))
  }
  mean.omega <- mean(omega)
  fcast <- return(list(t=tf,forecast=fcast, error=std.error))
}
