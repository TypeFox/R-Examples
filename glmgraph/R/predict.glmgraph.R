predict.glmgraph <- function(object, X,type=c("response", "coefficients", "class", "nzeros","link"),lambda1,lambda2,...) {
  type <- match.arg(type)
  betas <- coef(object,lambda1, lambda2,...)
  
  if (type=="coefficients") return(betas)
  if (type=="nzeros") {
  		tmp <- lapply(betas,function(beta) apply(beta[-1,,drop=FALSE]!=0,2,sum))
  		names(tmp) <- names(betas)
  		return(tmp)
  }
  eta <- mu <-  list()
  nlambda2 <- length(betas)
  
  for(i in 1:nlambda2) eta[[i]] <- sweep(X %*% betas[[i]][-1,,drop=FALSE], 2, betas[[i]][1,,drop=FALSE], "+")
  names(eta) <- names(betas)
  
  if(type=="link") return(eta);
  if (type=="response") {
  	 if (object$family=="gaussian") return(eta)
	 if (object$family=="binomial") {
  	 	for(i in 1:nlambda2) mu[[i]] <- exp(eta[[i]])/(1+exp(eta[[i]]))
  	 	names(mu) <- names(betas)
  	 	return(mu)
  	 }
  }
  
  if (type=="class") {
    if (object$family=="binomial") {
    	tmp <- lapply(eta,function(x) (x>0)*1)
    	names(tmp) <- names(betas)
      	return(tmp)
    } else {
      stop("type='class' can only be used with family='binomial'")
    }
  }
}








