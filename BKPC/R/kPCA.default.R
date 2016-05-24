kPCA.default <-
function(x, theta = NULL, ...) {
  
  if (is.null(x)) stop("error: no data supplied!")
  gaussKernel <- gaussKern(x, theta = theta)  
  K <- gaussKernel$K
  theta <- gaussKernel$theta
  
  
  objectPC <- getPrincipalComponents(K)
  object <- list(KPCs = objectPC$KPCs, Es = objectPC$Es, Vecs = objectPC$Vecs, K = K, theta = theta, x  = x)
  class(object) = "kPCA"  
  return(object)
}
