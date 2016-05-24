predict.kPCA <-
function(object, newdata, ...){
  
  if (inherits(object,  "kPCA.kern")){
    if (inherits(object,  "kPCA.kernelMatrix")){
      if (inherits(newdata,  "kernelMatrix")) testK <- newdata
      else stop("error: newdata should be of class 'kernelMatrix'") 
    }
    else{
      if (inherits(newdata,  "kern")) testK <- newdata
      else stop("error: newdata should be of class 'kern'") 
    }
  }
  
  if (!inherits(object,  "kPCA.kern")){
    if (inherits(newdata,  "kern"))stop("error: newdata is of class 'kern'") 
    else if (inherits(newdata,  "kernelMatrix"))stop("error: newdata is of class 'kernelMatrix'") 
    else  testK <- gaussKern(object$x, newdata, object$theta)$K
  }
  
  nr <- dim(testK)[1]
  nc <- dim(testK)[2]
  A <- matrix(1/nc, nc, nc)
  A2 <- matrix(1/nc, nr, nc)
  testKtild <- testK - A2 %*% object$K - testK %*% A + A2 %*% object$K %*% A  
  projected <- testKtild %*% object$Vecs
  colnames(projected) <- colnames(projected, do.NULL = FALSE, prefix = "KPC.")
  
  return(projected)
}
