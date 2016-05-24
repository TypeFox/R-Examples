print.LMmixed <-function(x, ...){
cat("Call:\n")
print(x$call)
cat("\nConvergence info:\n")
print(cbind(LogLik=x$lk,np=x$np,BIC=x$bic))
}