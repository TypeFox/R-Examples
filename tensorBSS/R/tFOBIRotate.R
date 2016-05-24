tFOBIRotate <-
function(x){
  r <- length(dim(x)) - 1
  
  rotateStack <- list()
  for(m in 1:r){
    rotateStack[[m]] <- t(eigenVectors(mModeTFOBIMatrix(x, m)))
  }
  
  
  for(m in 1:r){
    x <- tensorTransform(x, rotateStack[[m]], m)
  }
  
  return(list(x = x, U = rotateStack))
}
