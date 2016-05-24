computeOutput1 <- function(x,weights) {
  
  H <- length(weights@alpha_h)
  
  # calculate interim vector r
  r <- c(1:H)
  r <- vapply(r, function(h) weights@alpha_h[h] + sum(weights@w_ih[ ,h]*x),1)
  
  # calculate interim value z
  z <- weights@alpha + sum(weights@w_h*logistic(r))
  
  return(logistic(z))
}

