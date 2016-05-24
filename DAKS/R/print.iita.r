print.iita<-function(x, ...){
cat("\n \t Inductive Item Tree Analysis\n")
cat("\nAlgorithm:")
if(x$v == 1){cat(" minimized corrected IITA\n")}
if(x$v == 2){cat(" corrected IITA\n")}
if(x$v == 3){cat(" original IITA\n")}
cat("\nquasi order: ")
print(x$implications)
}
