LogPq <- function(x) {
  #For x = logit(p), this function returns log(p) and log(1-p). 
  #Special care is taken to ensure accuracy when coordinates 
  #are numerically close to 0 or 1.
  
  #Check inputs
  if (!exists("x")) {
    stop("x is not defined")
  }
  if (!is.numeric(x)) {
    stop("x is not numeric")
  }
  
  #Condition to choose which method is most numerically accurate
  ok <- (x < 0)
  
  #Make empty vectors
  logp <- logq <- rep(NA, length(x))
  
  #Find p/(1 - p)
  ex <- exp(x)
  
  #Calculate logp and logq for x < 0
  if (any(ok)) { 
    lq <- -log1p(ex[ok])
    logp[ok] <- lq + x[ok]
    logq[ok] <- lq
  }
  
  #Calculate logp and logq for x >= 0
  if (any(!ok)) {
    lp <- -log1p(1/ex[!ok])
    logp[!ok] <- lp
    logq[!ok] <- lp - x[!ok]
  }
  
  #Return logp and logq
  return(list(logp = logp, logq = logq))
}
