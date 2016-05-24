solvePR <-
function(mean, variance){
  p <- 1 - (mean/variance)
  r <- (mean)^2/(variance-mean)
  return(c(p,r))
}
