summary.iita<-function(object, ...){
cat("\n \t Inductive Item Tree Analysis\n")
cat("\nAlgorithm:")
if(object$v == 1){cat(" minimized corrected IITA\n")}
if(object$v == 2){cat(" corrected IITA\n")}
if(object$v == 3){cat(" original IITA\n")}
cat("diff values:",round(object$diff, digits = 3), " ")
cat("\nquasi order: ")
print(object$implications)
cat("error rate: ")
write(round(object$error.rate, digits= 3 ), "")
cat("index in the selection set: ")	
write(object$selection.set.index, "")
}
