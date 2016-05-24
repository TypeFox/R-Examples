summary.est_multi_poly_within <-function(object, ...){
	
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
    cat("\nClass weights for the 1st latent variable:\n")
    print(round(out$piv1,4))
    cat("\nClass weights for the 2nd latent variable:\n")
    print(round(out$piv2,4))
    cat("\nConditional response probabilities:\n")
    print(round(out$Phi,4))
    cat("\nEstimated abilities for the 1st latent variable:\n")
    print(round(out$Th1,4))    
    cat("\nEstimated abilities for the 2nd latent variable:\n")
    print(round(out$Th2,4))            
    cat("\n")

}