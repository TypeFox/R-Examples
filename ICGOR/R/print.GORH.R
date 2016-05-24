print.GORH <-
function(x,...){
   cat("Call:\n")
   print(x$call)
   cat("\nCoefficients:\n")
   tem.x<-x$ParEst$Beta
   names(tem.x)<-all.vars(x$survfun)[-c(1:2)]
   print(tem.x)
   cat("\nCovariance Matrix:\n")
   tem.v<-x$ParVcov$vcov.bg
   colnames(tem.v)<-rownames(tem.v)<-all.vars(x$survfun)[-c(1:2)]
   print(tem.v)
   print(paste("Loglik=",round(x$ParEst$loglik,2)))
   class(x)<-"GORH"
}
