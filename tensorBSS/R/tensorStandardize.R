tensorStandardize <-
function(x){
  r <- length(dim(x)) - 1
  
  x <- tensorCentering(x)
  invCovStack <- list()
  for(m in 1:r){
    invCovStack[[m]] <- symmetricPower(mModeCovariance(x, m), -0.5)
  }
  
  for(m in 1:r){
    x <- tensorTransform(x, invCovStack[[m]], m)
  }
  
  return(list(x = x, S = invCovStack))
}
