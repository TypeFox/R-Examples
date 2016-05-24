varComprob.fit <- function(Y, X, V, control=varComprob.control(), ...) {
# Y: matrix. dim(y)=c(p,n) 
# X: array. dim(x)=c(p,n,k)
# V: array. dim(V)=c(p,p,R)  
  init <- varComprob.initial(y=Y, x=X, V=V, beta=control$init$beta, gamma=control$init$gamma, eta0=control$init$eta0, Sigma=control$init$Sigma, scales=control$init$scales, scale=control$init$scale, control=control, ...)
  if (control$method=="compositeS")
    ans <- varComprob.compositeS(y=Y, x=X, V=V, beta=init$beta, gamma=init$gamma, eta0=init$eta0, scales=init$scales, control=control, ...)
  if (control$method=="compositeTau")
    ans <- varComprob.compositeTau(y=Y, x=X, V=V, beta=init$beta, gamma=init$gamma, eta0=init$eta0, scales=init$scales, control=control, ...)
  if (control$method=="S")
    ans <- varComprob.S(y=Y, x=X, V=V, beta=init$beta, gamma=init$gamma, eta0=init$eta0, scale=init$scale, control=control, ...)
  class(ans) <- c(class(ans), 'varComprob.fit')
  return(ans)
}
