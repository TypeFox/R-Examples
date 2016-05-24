QQdiff=function(object,items,plot=2, breaks=15,quant=NULL,maxRT=NULL){
  nit=object$nit
  N=object$N
  ai=object$par.log[1:nit]
  vi=object$par.log[(nit+1):(2*nit)]
  ter=object$par.log[(2*nit+1):(3*nit)]
  sd2A=object$par.log[3*nit+1]
  sd2V=object$par.log[3*nit+2]
  model=object$model
  nq=object$nq
  W=object$W
  A=object$A
  x=object$score
  rt=object$rt
  eps=object$control$eps
   
  if(ncol(x)!=nit | nrow(x)!=N) stop("The data and the diffIRT object do not match.\n")
  if(object$conv!=0) stop("The solution in 'object' has not converged.\n")
  if(object$nr_fail!=0) cat("Warning: The solution provided in 'object' might be instable.\n")
  if(plot!=1 & plot!=2) stop("Argument 'plot' should be either 1 or 2.\n")
  if(sum((items>object$nit)*1)!=0) stop("Non existent item number(s) provided.\n")
  if(!is.null(quant)) if(quant>N) stop("Number of quantiles exceeds the number of subjects.\n")
  if(!is.null(maxRT)) if(length(maxRT)!=length(items)) 
    stop("If you provide 'maxRT', it should be of the same size as specified in the 'items' argument") 
  if(plot==2) par(mfrow=c(ceiling(length(items)),2))
  if(plot==1) par(mfrow=c(ceiling(length(items)/2),2))
  
  if(is.null(quant)) qq=N else qq=quant
  prob=seq(0.05,.95,length=qq)
  
  outO=matrix(,qq,length(items))
  outE=matrix(,qq,length(items))
  for(i in 1:length(items)){
    if(plot==2){
      x=seq(ter[items[i]],max(rt[,items[i]],na.rm=T),length=1000)
      pred=mdiffIRT(x,ai[items[i]],vi[items[i]],ter[items[i]],sd2A,sd2V,A,W,model,eps)
      hist(rt[,items[i]],freq=F,main=paste("item",items[i]),xlab="observed and predicted distribution",ylim=c(0,max(pred)),breaks)
      lines(x,pred,lty=2)
    }
    outO[,i]=quantile(rt[,items[i]],probs=prob,na.rm=TRUE)
    if(!is.null(maxRT)) mRT=maxRT[i]
    else mRT=2*max(rt[,items[i]],na.rm=T)
    if(is.na(mRT)) mRT=4*max(rt[,items[i]],na.rm=T)
    outE[,i]=qdiffIRT(prob,ai[items[i]],vi[items[i]],ter[items[i]],sd2A,sd2V,maxRT=mRT,A,W,model,eps)
    plot(outE[,i],outO[,i],xlab="theoretical quantity",ylab="observed quantitiy",main=paste("QQ plot item",items[i]))
    abline(0,1)
  }
  res=list(qexp=outE,qobs=outO)
  invisible(res)

}



