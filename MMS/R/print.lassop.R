print.lassop <-
function(x,...) {
    # 0. part
    print(x$call)

    # 2a. part
    cat("","\n")
    for(i in 1:length(x$mu))
{   cat("Results for mu=",x$mu[i],"\n")
    cat("Number of relevant fixed-effects:",sum(x$beta!=0),"\n")
  	cat("which are:",which(x$beta!=0),"\n\n")
  	cat("\n")
  	
cat("Variance of the random effects:\n")
print(x$Psi)
cat("\n")
cat("Variance of the residuals", x$sigma_e,"\n")
#cat("-------------------------------\n\n")
}
}



