print.GORMC <-
function(x,...){
   cat("Call:\n")
   print(x$call)
   cat("\nCoefficients:\n")
   cat("\nIncidence Part:\n")
   tem.z<-x$ParEst$Eta
   names(tem.z)<-c("Intercept",all.vars(x$curefun))
   print(tem.z)
   cat("\nLatency Part:\n")
   tem.x<-x$ParEst$Beta
   names(tem.x)<-all.vars(x$survfun)[-c(1:2)]
   print(tem.x)
   cat("\nCovariance Matrix:\n")
   tem.v<-x$ParVcov$vcov.bg
   colnames(tem.v)<-rownames(tem.v)<-c(paste("INC:",c("Intercept",all.vars(x$curefun))),paste("LAT:",all.vars(x$survfun)[-c(1:2)]))
   print(tem.v)
   print(paste("Loglik=",round(x$ParEst$loglik,2)))
   class(x)<-"GORMC"
}
