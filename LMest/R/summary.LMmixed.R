summary.LMmixed<-function(object,...){ 
cat("Call:\n")
print(object$call)
out_se = object$call$out_se
cat("\nCoefficients:\n")
cat("\nMass probabilities:\n")
print(round(object$la,4)) 
if(!is.null(object$sela)){
	cat("\nStandard errors for the mass probabilities:\n")
	print(round(object$sela,4))
} 
cat("\nInitial probabilities:\n")
print(round(object$Piv,4))
if(!is.null(object$sePiv)){
	cat("\nStandard errors for the initial probabilities:\n")
	print(round(object$sePiv,4))
}
cat("\nTransition probabilities:\n")
print(round(object$Pi,4))
if(!is.null(object$sePi)){
	cat("\nStandard errors for the transition probabilities:\n")
	print(round(object$sePi,4))
}
cat("\nConditional response probabilities:\n")
print(round(object$Psi,4))
if(!is.null(object$sePsi)){
	cat("\nStandard errors for the conditional response probabilities:\n")
	print(round(object$sePsi,4))
}
}