tSOBIRotate <-
function(x, lags, maxiter, eps){
  r <- length(dim(x)) - 1
  
  rotateStack <- list()
  for(m in 1:r){
    pm <- dim(x)[m]
    lagStack <- matrix(0, pm*length(lags), pm)
    for(i in 1:length(lags)){
        lagmat <- mModeAutoCovariance(x, m, lags[i], center=FALSE)
        lagStack[((i-1)*pm + 1):(i*pm) , 1:pm] <- (lagmat + t(lagmat))/2
    }
    rotateStack[[m]] <- t(frjd(lagStack, maxiter=maxiter, eps=eps)$V)
  }
  
  
  for(m in 1:r){
    x <- tensorTransform(x, rotateStack[[m]], m)
  }
  
  return(list(x = x, U = rotateStack))
}
