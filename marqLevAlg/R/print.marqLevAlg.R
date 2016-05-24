print.marqLevAlg <- function(x,digits=8,...){
if (!inherits(x, "marqLevAlg")) stop("use only with \"marqLevAlg\" objects")

cl <- x$cl
cat(" \n")
dput(cl)
cat(" \n")
cat("                   Robust marqLevAlg algorithm                   ", "\n")
cat(" \n")
cat("Number of parameters:", length(x$b)," \n")
cat(" \n")
cat("Iteration process:", "\n")

if(x$istop==1) cat("      Convergence criteria satisfied","\n")
if(x$istop==2) cat("      Maximum number of iteration reached without convergence","\n")
if(x$istop==4|x$istop==5)  {
	cat("      The program stopped abnormally. No results can be displayed.\n")
}else{
cat(" \n")
cat("Values:", "\n")
cat(" \n")
id <- 1:length(x$b)
indice <- rep(id*(id+1)/2)
se <-sqrt(x$v[indice])
wald <- (x$b/se)**2
z <- abs(qnorm((1 + .95)/2))

tmp <- data.frame("coef"=format(round(x$b,digits)),"SE coef"=format(round(se,digits)),"Wald"=format(wald,4),"P-value"=format.pval(1 - pchisq(wald, 1),digits=digits,eps=0.0001))
print(tmp,row.names=F)
cat(" \n")
cat("Number of iterations: ", x$ni, "\n")
cat(" \n")
cat("Convergence criteria: parameters stability=", round(x$ca[1],digits), "\n")
cat("                    : likelihood stability=", round(x$cb,digits), "\n") 
if (x$ier == -1){
	cat("                    : Matrix inversion for RDM failed \n")	
}else{
	cat("                    : Matrix inversion for RDM successful \n")
}
cat("                    : relative distance to maximum(RDM)=", round(x$rdm,digits), "\n")
cat(" \n")
cat("Goodness-of-fit statistics:", "\n")
cat("      minimum log-likelihood:", round(x$fn.value,digits)," \n")
cat(" \n")
cat(" \n")
cat(" \n")

}
}
