`leaveOneOutFun` <-
function(param, model, envir=NULL) {
	    
	if (model@case == "LLconcentration_beta_sigma2") {
		
		model@covariance <- vect2covparam(model@covariance, param)
		model@covariance@sd2 <- 1		# to get the correlation matrix
		
	  R <- covMatrix(model@covariance, model@X)[[1]]
		T <- chol(R)
		    
   	M <- backsolve(t(T), model@F, upper.tri = FALSE)
		
    Rinv <- chol2inv(T)             # cost : n*n*n/3
		if (model@known.param=="None"){
      Rinv.F <- Rinv %*% (model@F)    # cost : 2*n*n*p
		  T.M <- chol(crossprod(M))       # cost : p*p*p/3, neglected
		  aux <- backsolve(t(T.M), t(Rinv.F), upper.tri=FALSE)   # cost : p*p*n, neglected
		  Q <- Rinv - crossprod(aux)      # cost : 2*n*n*(p-1/2)
      Q.y <- Q %*% (model@y)          # cost : 2*n*n
		} else if (model@known.param=="Trend") {
      Q <- Rinv
      Q.y <- Q %*% (model@y - model@F %*% model@trend.coef)
		}
		sigma2LOO <- 1/diag(Q)
		errorsLOO <- sigma2LOO * (Q.y)       # cost : n, neglected 
		
    LOOfun <- as.numeric(crossprod(errorsLOO)/model@n)
    
		if (!is.null(envir)) { 
      envir$R <- R
      envir$Q <- Q
      envir$Q.y <- Q.y
      envir$errorsLOO <- errorsLOO
      envir$sigma2LOO <- sigma2LOO
		}
		
	} else {
		stop("leave-One-Out is not available for this model")
	} 
	
	return(LOOfun)
}