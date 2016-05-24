bkpc.kernelMatrix <-
function(x, y, n.kpc = NULL, thin = 100, 
                              n.iter = 100000, std = 10, g1 = 0.001, 
                              g2 = 0.001, g3 = 1, g4 = 1, initSigmasq = NULL, initBeta = NULL, initTau = NULL, intercept = TRUE, rotate = TRUE, ...){
  
  old.class <- class(x)
  class(x) <- "kern"
  object <- bkpc.kern(x, y, n.kpc, thin, n.iter, std, g1, g2, g3, g4, initSigmasq, initBeta, initTau,  intercept, rotate)
  object$x <- x
  class(object$x) <- old.class
  
  if(rotate){
   class(object$kPCA$K) <- old.class
   class(object$kPCA) <- c("kPCA.kernelMatrix", class(object$kPCA)) 
  }
  class(object) = c("bkpc.kernelMatrix", class(object)) 
  return(object)
}
