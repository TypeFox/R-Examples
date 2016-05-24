print.Bayesian.SGoF <-
function(x, ...){
cat("Call:\n") 
print(x$call) 
cat("\n")
cat("Parameters:","\n")
cat("alpha=",x$alpha,"\n")
cat("gamma=",x$gamma,"\n")
cat("P0=",x$P0,"\n")
cat("a0=",x$a0,"\n")
cat("b0=",x$b0,"\n")
cat("\n")
cat("Rejections:\n") 
print(x$Rejections)
}
