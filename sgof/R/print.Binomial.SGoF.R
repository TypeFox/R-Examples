print.Binomial.SGoF <-
function(x, ...){
cat("Call:\n") 
print(x$call) 
cat("\n")
cat("Parameters:","\n")
cat("alpha=",x$alpha,"\n")
cat("gamma=",x$gamma,"\n")
cat("\n")
cat("Rejections:\n") 
print(x$Rejections)
}
