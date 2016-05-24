print.Discrete.SGoF <-
function(x,...){
cat("Call:\n") 
print(x$call) 
cat("\n")
cat("Parameters:","\n")
cat("alpha=",x$alpha,"\n")
cat("gamma=",x$gamma,"\n")
if(is.na(x$K)==FALSE){cat("K=",x$K,"\n")}
if(is.na(x$method)==FALSE){cat("Method=",x$method,"\n")}
cat("Discrete=",x$Discrete,"\n")
cat("Sides=",x$Sides,"\n")
cat("\n")
cat("Rejections:\n") 
print(x$Rejections)
cat("FDR:\n") 
print(x$FDR)
}
