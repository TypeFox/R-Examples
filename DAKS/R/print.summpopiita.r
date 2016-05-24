print.summpopiita<-function(x, ...){
cat("\n \t Inductive Item Tree Analysis in population values\n")
cat("\nAlgorithm:")
if(x$v == 1){cat(" minimized corrected IITA\n")}
if(x$v == 2){cat(" corrected IITA\n")}		
if(x$v == 3){cat(" original IITA\n")}		
cat("\npopulation diff values:\n")
print(round(x$pop.diff, digits = 3))
cat("\npopulation error rates:\n")
print(round(x$error.pop, digits = 3))
cat("\npopulation matrix:\n")
print(round(x$pop.matrix, digits =3))
cat("\nobtained selection set:\n")
print(x$selection.set)
}
