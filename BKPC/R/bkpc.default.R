bkpc.default <-
function(x, y,  theta = NULL, n.kpc = NULL, thin = 100, 
                         n.iter = 100000, std = 10,  g1 = 0.001, 
                         g2 = 0.001,  g3 = 1, g4 = 1, initSigmasq = NULL, initBeta = NULL, initTau = NULL,  intercept = TRUE, rotate = TRUE, ...) {


  gaussianKernel <-  gaussKern(x, theta = theta) 
  K <- gaussianKernel$K # class kern  
    
  object <- bkpc.kern(K, y, n.kpc, thin, n.iter, std, g1, g2, g3, g4, initSigmasq, initBeta, initTau, intercept, rotate)
  
  object$x <- x
  object$theta <- gaussianKernel$theta  
  
  class(object) <- "bkpc"
  return(object)
  
}
