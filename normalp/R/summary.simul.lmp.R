summary.simul.lmp<-function(object,...){
ans<-object
szca<-nrow(ans$dat)
ans$szca<-szca
class(ans)<-"summary.simul.lmp"
ans
}
  
print.summary.simul.lmp<-function(x,...){
a<-x
cat("Results:\n")
print(a$table)
cat("\nCoefficients:\n")
print(a$par)
cat("\nFormula:\n")
print(a$frm)
cat("\nNumber of samples:\n")
print(as.symbol(a$szca))
cat("\nValue of p:\n")
print(as.symbol(a$p))
if(x$lp==FALSE){
cat("\nNumber of samples with problems on convergence\n")
print(as.symbol(a$iter))
}
cat("\n")
invisible(a)
}

