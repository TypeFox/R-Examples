summary.test_dim <-function(object, ...){
	
# preliminaries
	out0 = object$out0
	out1 = object$out1
	out = object
# print output	
	cat("\nCall:\n")
    print(out$call)
    cat("\nTesting dimension output:\n")
    table = c(round(out0$lk,3),round(out0$aic,3),round(out0$bic,3),round(out0$np,0),
              round(out1$lk,3),round(out1$aic,3),round(out1$bic,3),round(out1$np,0),
              round(out$dev,3),out$df,round(out$pv,3))
    table = matrix(table,11,1)
    colnames(table) = ""
    rownames(table) = c("Log-likelihood of the constrained model","AIC of the constrained model",
                        "BIC of the constrained model","N.parameters of the constrained model",
                        "Log-likelihood of the unconstrained model","AIC of the unconstrained model",
                        "BIC of the unconstrained model","N.parameters of the unconstrained model",
                        "Deviance","Degrees of freedom","p-value")
	print(table)	
    cat("\n")

}