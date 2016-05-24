summary.diffIRT = function(object,digits=3,...){
  object->x
  summ=matrix(x$par[1:(x$nit*3)],x$nit,3)
  if(x$se){ 
    se=matrix(x$std.err[1:(x$nit*3)],x$nit,3)
    summ=cbind(summ[,1],se[,1],summ[,2],se[,2],summ[,3],se[,3])
    dimnames(summ)[[2]]<-c("a[i]","se.a[i]","v[i]","se.v[i]","Ter[i]","se.Ter[i]")
  } else dimnames(summ)[[2]]<-c("a[i]","v[i]","Ter[i]")
  dimnames(summ)[[1]]<-paste("item ",1:x$nit,sep="")
  summ=round(summ,digits,...)
  cat("\n        RESULTS\n")
  cat("        -------\n")
  print(summ)
  sdA=x$par[(x$nit*3+1)]
  sdV=x$par[(x$nit*3+2)]
  VarAp=(exp(sdA^2)-1)*exp(sdA^2)
  if(x$model==1) VarVp=sdV^2
  if(x$model==2) VarVp=(exp(sdV^2)-1)*exp(sdV^2)
  cat("\n")
  cat("omega[gamma]:",round(sdA,digits),"\n")
  if(x$se) cat("std. err:    ", round(x$std.err[x$nit*3+1],digits),"\n")
  cat("omega[theta]:",round(sdV,digits),"\n")
  if(x$se) cat("std. err:    ",round(x$std.err[x$nit*3+2],digits),"\n")
  if(x$se) if((TRUE %in% (x$std.err==-999))) cat("\n* -999: this parameter is fixed\n") else  cat("\n")
  cat("---------------------------\n")
  cat("\n")
  cat("--POPULATION DESCRIPTIVES--\n")
  cat("Person Boundary (lognormal): Var(gamma): ",round((VarAp),digits),"\n")
  if(x$model==1) cat("Person Drift (normal): Var(theta): ",round((VarVp),digits),"\n")
  if(x$model==2) cat("Person Drift (lognormal):    Var(theta): ",round((VarVp),digits),"\n")
  cat("---------------------------\n")
  cat("\n")
  cat("---MODEL FIT STATISTICS---\n")
  cat("-2 x logLikelihood:", round(x$totLL,3), "\n")
  cat("\nno. of parameters:", x$npar,"\n")
  cat("AIC: ",round(x$AIC,3),"\n")
  cat("BIC: ",round(x$BIC,3),"\n")
  cat("sBIC:",round(x$sBIC,3),"\n")
  cat("DIC: ",round(x$DIC,3),"\n")
  cat("---------------------------\n")
 }
 