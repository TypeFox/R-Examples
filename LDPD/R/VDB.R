VDBCAPExpAp <- function(k, x) {
  # One parameter approximation of cAP curve
  # Args:
  #   k:      approximation parameters (describes convexity of CAP curve)
  #   x:      point of borrower's unconditional cumulative distribution
  # Returns:
  #           value of approximated CAP curve in a point x
  return ((1 - exp(-k * x)) / (1 - exp(-k)))
}

VDBCAPDer <- function(x, k) {
  # One parameter approximation of Derivative ofcAP curve
  # Args:
  #   k:      approximation parameters (describes convexity of CAP curve)
  #   x:      point of borrower's unconditional cumulative distribution
  # Returns:
  #                   value of approximated Derivative of CAP curve in a point x
  return  ((k * exp(-k * x)/(1 - exp(-k))))
}

VDBAR <- function(k, pd.uncond) {
  # Value of AR given approximation parameter and unconditional PD of the sample
  # Args:
  #   k:          approximation parameters (describes convexity of CAP curve)
  #   pd.uncond:  unconditional PD of the sample 
  # Returns:
  #               Value of AR given approximation parameter and unconditional PD of the sample
  return((2 * (1 / (1 - exp(-k)) - 1 / k - 0.5) / (1 - pd.uncond)))
}

VDBGetK <- function(AR, pd.uncond) {
  # Finds an approximation parameter k for CAP curve
  #   k:          approximation parameters (describes convexity of CAP curve)
  #   pd.uncond:  unconditional PD of the sample 
  # Returns:
  #               approximation parameters (describes convexity of CAP curve)
  k <- uniroot(function(x) {VDBAR(x, pd.uncond) - AR}, lower = 1e-5, upper = 100)$root
  return(k)
}
VDBCalibratePD <- function(portf.uncond, pd.uncond.old, pd.uncond.new, AR, rating.type) {
  # Calibrates PDs acording to approach proposed by M. van der Burgt.
  # Args:
  #   portf.uncond:     unconditional portfolio distribution from the worst to the best credit quality
  #   pd.uncond.new:    target Mean PD (Central Tendency) for the porfolio  
  #   pd.uncond.old:    unconditional PD of the sample on wich AR had been estimated
  #   AR:               AR of the ranking model
  #   rating.type:      In case RATING, each item in the portf.uncond contains number of companies in a given rating class
  # Returns:
  #   lambda:           convexity parameter of the calibration curve
  #   pd.cond:          conditional PDs after calibration 
  #   portf.cumdist:    cumulative portfolio distribution needed to estimate logit PDs (conditional on non-default if such data is given)
  #   portf.uncond:     unconditional portfolio distribution from the worst to the best credit quality
  #   rating.type:      In case RATING, each item in the portf.uncond contains number of companies in a given rating class  
  
  if (rating.type == 'RATING') { 
    portf.cum <- cumsum(portf.uncond)
    portf.dist <- (portf.cum + c(0, portf.cum[-length(portf.cum)])) / (2 * sum(portf.uncond))
  } else {
    portf.dist <- ecdf(portf.uncond)(portf.uncond)
    portf.dist[length(portf.dist)] <- (1 + portf.dist[length(portf.dist) - 1]) / 2 
  }
  k <- VDBGetK(AR, pd.uncond.old)

  pd.cond <- pd.uncond.new * sapply(portf.dist, VDBCAPDer, k = k)
  
  rez <- list()
  rez$lambda            <- k  
  rez$condPD           <- pd.cond
  rez$portf.cumdist    <- portf.dist
  rez$portf.uncond     <- ifelse(rating.type == 'RATING', 
                                 t(portf.uncond / sum(portf.uncond)) %*% pd.cond, 
                                 mean(pd.cond))
  rez$rating.type      <- rating.type
  return(rez)
  
}


