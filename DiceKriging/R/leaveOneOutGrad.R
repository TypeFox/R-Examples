`leaveOneOutGrad` <-
function(param, model, envir) {
	
  if (model@case == "LLconcentration_beta_sigma2") {
    
		R <- envir$R
    Q <- envir$Q
    Q.y <- envir$Q.y
    errorsLOO <- envir$errorsLOO
    sigma2LOO <- envir$sigma2LOO

		model@covariance <- vect2covparam(model@covariance, param)
		model@covariance@sd2 <- 1		# to get the correlation matrix
	
		nparam <- length(param)
		LOOfunDer <- matrix(0, nparam, 1)
										
		for (k in 1:nparam) {
			gradR.k <- covMatrixDerivative(model@covariance, X=model@X, C0=R, k=k)
			diagdQ <- - diagABA(A=Q, B=gradR.k)
			dsigma2LOO <- - (sigma2LOO^2) * diagdQ
			derrorsLOO <- dsigma2LOO * Q.y - sigma2LOO * (Q%*%(gradR.k%*%Q.y))
			LOOfunDer[k] <- 2*crossprod(errorsLOO, derrorsLOO)/model@n
		}
	
	} else {
    stop("leave-One-Out is not available for this model")
	} 
	
	return(LOOfunDer)
}