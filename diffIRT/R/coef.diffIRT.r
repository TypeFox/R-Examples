coef.diffIRT = function(object,...){ 
  summ=matrix(object$par[1:(object$nit*3)],object$nit,3)
  dimnames(summ)[[2]]<-c("ai","vi","ter")
  dimnames(summ)[[1]]<-paste("item ",1:object$nit,sep="")
  sdA=object$par[(object$nit*3+1)]
  sdV=object$par[(object$nit*3+2)]
  pop=matrix(c(sdA,sdV),1,2)
  dimnames(pop)[[2]]=c("omega[gamma]","omega[theta]")
  return(list(item=summ,pop=pop))
}


