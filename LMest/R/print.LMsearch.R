print.LMsearch <-function(x,...){
cat("Call:\n")
print(x$call)
#print(rbind(lkv=x$lkv,aicv=x$aicv,bicv=x$bicv))
cat("lk = ",x$lkv,"\n")
cat("aic = ",x$aicv,"\n")
cat("bic = ",x$bicv,"\n")
}