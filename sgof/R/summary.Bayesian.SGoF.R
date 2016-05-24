summary.Bayesian.SGoF <-
function(object, ...){
cat("Call:\n")
 print(object$call)
cat("\n")
cat("Parameters:","\n")
cat("alpha=",object$alpha,"\n")
cat("gamma=",object$gamma,"\n")
cat("P0=",object$P0,"\n")
cat("a0=",object$a0,"\n")
cat("b0=",object$b0,"\n")

cat("\n")
res <- list(Rejections=object$Rejections,FDR=object$FDR,Posterior=object$Posterior,s=object$s,s.alpha=object$s.alpha )


class(res) <- "summary.Bayesian.SGoF"
return(res)
}
