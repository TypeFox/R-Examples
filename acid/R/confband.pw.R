confband.pw <-
function(samples, level=0.95){
  # at least k functions completely cotained in the band
  n <- nrow(samples)
  p <- ncol(samples)
  
  k <- trunc(level*n)+1
  
  pmean <- apply(samples, 2, mean)
  pqlow <- apply(samples, 2, quantile, probs=(1-level)/2) #quantile for each column, i.e. each point in line
  pqupp <- apply(samples, 2, quantile, probs=1-(1-level)/2)
  lower <- pqlow
  upper <- pqupp
  
  return(list(lower=lower, upper=upper))
}
