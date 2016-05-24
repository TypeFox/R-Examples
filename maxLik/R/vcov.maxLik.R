## maxLik
vcov.maxLik <- function(object, eigentol=1e-12, ...) {
   ## if exists $varcovar, take it
   if(!is.null(object$varcovar))
       return(object$varcovar)
   ## otherwise invert hessian
   activePar <- activePar(object)
   if(!is.null(hess <- hessian(object))) {
      hess <- hessian(object)[activePar, activePar,drop=FALSE] 
      hessev <- abs(eigen(hess, symmetric=TRUE, only.values=TRUE)$values)
      varcovar <- matrix(0, nParam.maxim(object), nParam.maxim(object))
                           # this makes the fixed parameters to 0
      rownames( varcovar ) <- colnames(varcovar ) <- names(coef.maxLik(object))
      if(min(hessev) > (eigentol*max(hessev))) {
      ## If hessian is not singular, fill in the free parameter values
         varcovar[activePar,activePar] <- solve(-hessian(object)[activePar,activePar])
         # guarantee that the returned variance covariance matrix is symmetric
         varcovar <- ( varcovar + t( varcovar ) ) / 2
      }
      else {
      ## If singular, the free parameter values will be Inf
         varcovar[activePar,activePar] <- Inf
      }
      return(varcovar)
   }
   else
       return(NULL)
}
