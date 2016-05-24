summary.est_multi_poly <-function(object, ...){
	
# preliminaries
	out = object
# print output	
	cat("\nCall:\n")
    print(out$call)
    cat("\nLog-likelihood:\n")
    print(round(out$lk,2))
    cat("\nAIC:\n")
    print(round(out$aic,2))
    cat("\nBIC:\n")
    print(round(out$bic,2))
    cat("\nClass weights:\n")
    print(round(out$piv,4))
    cat("\nConditional response probabilities:\n")
    print(round(out$Phi,4))    
    cat("\n")

}