summary.LMmanifest<-function(object,...){ 
	cat("Call:\n")
	print(object$call)
	cat("\nCoefficients:\n")
	cat("\n  Vector of cut-points:\n")
	print(round(object$mu,digits=4))
	cat("\n Support points for the latent states:\n")
	print(round(object$al,digits=4))
	cat("\n Estimate of the vector of regression parameters:\n")
	print(round(object$be,digits=4))	
	cat("\n Vector of initial probabilities:\n")
	print(round(object$la,digits=4))
	cat("\n Transition matrix:\n")
	print(round(object$PI,digits=4))
	
	if(is.null(object$si)==FALSE){
		cat("\n Sigma of the AR process:\n")
		print(round(object$si,digits=4))		
		cat("\n Parameter vector for AR:\n")
		print(round(object$rho,digits=4))
	}
	if(is.null(object$sebe)==FALSE){
		cat("\n Standard errors for the regression parameters:\n")
		print(round(object$sebe,digits=4))
	}
	if(is.null(object$selrho)==FALSE){
		cat("\n Standard errors for the logit type transformation of rho:\n")
		print(round(object$selrho,digits=4))
	}

	
}