########summary the class "Networks.STD"
summary.Networks.STD=function(object, ...){
  classification=sapply(1:ncol(object$trace),function(kk) return(mean(object$trace[,kk])))
  eids1=which(classification>0)
  eids0=which(classification<=0)
  classification[eids1]<-1
  classification[eids0]<-0
  cat("i.  Convergence results:","\n")
  print(object$convergence)
  cat("\n")
  cat("ii.  Hyper-Parameter Selection:","\n")
  cat("pi0=",unlist(object$model$parameter$pi0),"rho0=",unlist(object$model$parameter$rho0),"rho1=",unlist(object$model$parameter$rho1),"\n")
  cat("\n")
  cat("iii.  Classification Table:","\n")
  print(table(classification))
}
