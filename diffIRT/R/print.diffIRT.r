print.diffIRT = function(x,digits=3,...){
  x->object
  summ=matrix(object$par[1:(object$nit*3)],object$nit,3)
  dimnames(summ)[[2]]<-c("a[i]","v[i]","Ter[i]")
  dimnames(summ)[[1]]<-paste("item ",1:object$nit,sep="")
  summ=round(summ,digits)
  if(object$model==1) cat("\nRESULTS D-DIFFUSION IRT ANALYSES\n")
  if(object$model==2) cat("\nRESULTS Q-DIFFUSION IRT ANALYSES\n")
  cat("\nItem parameter estimates\n")
  cat("        -------\n")
  print(summ)
  cat("\n")
  sdA=object$par[(object$nit*3+1)]
  sdV=object$par[(object$nit*3+2)]
  cat("Population parameter estimates\n")
  cat("        -------\n")
  cat("omega[gamma]:",round(sdA,digits),"\n")
  cat("omega[theta]:",round(sdV,digits),"\n")
  cat("---------------------------\n")
}


