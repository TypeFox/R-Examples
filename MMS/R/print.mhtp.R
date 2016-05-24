print.mhtp <-
function(x,...) {
    # 0. part
    print(x$call)

    # 2a. part
    cat("","\n")
    alpha=as.numeric(colnames(x$fitted.values))  	     
    for(i in 1:length(alpha))
{   cat("Results for alpha=",alpha[i],"\n")
    cat("Number of relevant fixed-effects:",sum(x$beta[,i]!=0),"\n")
  	cat("which are:",which(x$beta[,i]!=0),"\n\n")
  	
cat("Variance of the random effects:\n")
print(x$Psi[,,i])
cat("\n")
cat("Variance of the residuals", x$sigma_e[i],"\n\n")
cat("-------------------------------\n\n")
}
}



