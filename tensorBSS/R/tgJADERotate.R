tgJADERotate <-
function(x, lags, maxiter, eps){
  r <- length(dim(x)) - 1
  
  rotateStack <- list()
  for(m in 1:r){
    pm <- dim(x)[m]
    lagStack <- NULL
    for(l in 1:length(lags)){
      for(i in 1:pm){
        for(j in 1:pm){
          lagStack <- rbind(lagStack, mModeTGJADEMatrix(x, m, i, j, lags[l], pm))
        }
      }
    }
    rotateStack[[m]] <- t(frjd(lagStack, maxiter=maxiter, eps=eps)$V)
  }
  
  
  for(m in 1:r){
    x <- tensorTransform(x, rotateStack[[m]], m)
  }
  
  return(list(x = x, U = rotateStack))
}
