simID <- function(n=1, obs=60, idi=2, cv2=0.5, level=NULL){
  # Simulator for Intermittent Demand Series
  #
  # Inputs:
  #   n           Number of time series to be generated.
  #   obs         Number of observation of each series.
  #   idi         Average intermittent demand interval of each series.
  #   cv2         Squared coefficient of variation of the non-zero demands.
  #   level       Mean level of the non-zero demands. If NULL, then a random level in [10,100] is selected.
  #
  # Outputs:
  #   series       A data matrix containing all the generated series.
  #
  # Example:
  #   dataset <- t(simID(100,60,idi=1.15,cv2=0.3))
  #
  # Notes:
  # This simulator assumes that non-zero demand arrivals follow a bernoulli distribution 
  # and the non-zero demands a negative binomial distribution
  # Based on the paper: Petropoulos F., Makridakis S., Assimakopoulos V. & Nikolopoulos K. (2014) 
  # "'Horses for Courses' in demand forecasting", European Journal of Operational Research, Vol. 237, No. 1, pp. 152-163
  #
  # Fotios Petropoulos, 2014 <fotpetr@gmail.com>

  series <- array(NA, c(obs, n))
  
  for (tsi in 1:n){
    if (is.null(level)){
      m <- runif(1,9,99)
    } else {
      m <- level - 1
    }
    
    if (cv2!=0){
      p <- (m/(cv2*((m+1)^2)))
      r <- m*p/(1-p)
      x <- rbinom(obs,1,1/idi) * (rnbinom(obs, r, p)+1) # rbern(obs, 1/idi) replaced by binomial
    } else {
      x <- rbinom(obs,1,1/idi) * round(m+1) # rbern(obs, 1/idi) replaced by binomial
    }
    
    series[, tsi] <- x
  }
  
  # Convert series to data.frame
  colnames(series) <- paste("ts",1:n)
  series <- data.frame(series)
  
  return(series)
  
}