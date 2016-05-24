crost.decomp <- function(data,init=c("naive","mean")){
  # Croston decomposition
  #
  # Inputs:
  #   data        Intermittent demand time series.
  #   init        Initial values for intervals. This can be:
  #                 x       - Numerical value;
  #                 "naive" - Initial interval is the first interval 
  #                           from start of time series;
  #                 "mean"  - Initial interval is the mean of all 
  #                           in sample intervals.
  #
  # Outputs:
  #   demand      Non-zero demand vector
  #   interval    Intervals vector
  #
  # Example:
  #   crost.decomp(ts.data1)
  #
  #
  # Nikolaos Kourentzes, 2015 <nikolaos@kourentzes.com>
  
  # Defaults
  init <- init[1]
  
  if (class(data) == "data.frame"){
    if (ncol(data)>1){
      warning("Data frame with more than one columns. Using only first one.")
    }
    data <- data[[1]]
  }
    
  n <- length(data)
  
  # Check number of non-zero values - need to have at least two
  if (sum(data!=0)<1){
    stop("Need to have at least on nonzero demand to decompose the time series.")
  }
  
  # Croston decomposition
  nzd <- which(data != 0)               # Find location on non-zero demand
  k <- length(nzd)
  z <- data[nzd]                        # Demand
  x <- c(nzd[1],nzd[2:k]-nzd[1:(k-1)])  # Intervals
  
  # Initialise intervals
  if (!(is.numeric(init))){
    if (init=="mean"){
      x[1] <- mean(x)
    }
  } else {
    x[1] <- init
  }
  
  return(list(demand=z, interval=x))
  
}