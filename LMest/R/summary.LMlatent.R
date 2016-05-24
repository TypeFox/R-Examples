summary.LMlatent<-function(object,...){ 
cat("Call:\n")
print(object$call)
cat("\nCoefficients:\n")
if(object$call$k>1){
	cat("\n Be - Parameters affecting the logit for the initial probabilities:\n")
	print(round(object$Be,4))
	if(is.null(object$seBe)==FALSE){
		cat("\n Standard errors for Be:\n")
		print(round(object$seBe,4))
	}
	if(object$param=="multilogit"){
	cat("\n Ga - Parameters affecting the logit for the transition probabilities:\n")
		print(round(object$Ga,4))	
		if(is.null(object$seGa)==FALSE){
			cat("\n Standard errors for Ga:\n")
			print(round(object$seGa,4))		
		}	
	}else if(object$param=="difflogit"){
		cat("\n Ga0 - Intercept affecting the logit for the transition probabilities:\n")
		print(round(object$Ga[[1]],4))
		cat("\n Ga1 - Regression parameters affecting the logit for the transition probabilities:\n")
		print(round(object$Ga[[2]],4))
		
		if(is.null(object$seGa)==FALSE){
			cat("\n Standard errors for Ga0:\n")
			print(round(object$seGa[[1]],4))
			cat("\n Standard errors for Ga1:\n")
			print(round(object$seGa[[2]],4))		
		}
	}
	
}
	cat("\n Psi - Conditional response probabilities:\n")
	print(round(object$Psi,4))
	if(is.null(object$sePsi)==FALSE){
		cat("\n Standard errors for the conditional response probability matrix:\n")
		print(round(object$sePsi,4))
	}
}