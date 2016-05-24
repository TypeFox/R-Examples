computeOutput2 <- function(x,weights) {
  
  H <- length(weights@alpha_2h)
  M <- length(weights@alpha_1m)
  
  # calculate interim vector r
  r <- c(1:M)
  r <- vapply(r, function(m) weights@alpha_1m[m] + sum(weights@w_im[ ,m]*x),1)
  
  # calculate interim value z
  z <- c(1:H)
  z <- vapply(z, function(h) weights@alpha_2h[h] + sum(weights@q_mh[ ,h]*logistic(r)),1)
  
  # calculate interim value s
  s <- weights@alpha + sum(weights@w_h*logistic(z))
  
  return(logistic(s))
}