cv.sparsenet=function(x,y,weights,type.measure=c("mse","mae"),...,nfolds=10,foldid,trace.it=FALSE){
this.call=match.call()
  type.measure=match.arg(type.measure)
  N=nrow(x)
  if(missing(weights))weights=rep(1.0,N)else weights=as.double(weights)

###Fit the model once to get dimensions etc of output
  y=drop(y) # we dont like matrix responses unless we need them
  sparsenet.object=sparsenet(x,y,...)
  parms=sparsenet.object$parms
  nz=sapply(predict(sparsenet.object,type="nonzero"),function(x)sapply(x,length))
  dd=dim(nz)
  ngamma=dd[2]
  nlams=matrix(0,nfolds,dd[2])
  predmat=array(NA,c(N,dd))
  if(missing(foldid)) foldid=sample(rep(seq(nfolds),length=N)) else nfolds=max(foldid)
  if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=10 recommended")
   outlist=as.list(seq(nfolds))
###Now fit the nfold models and store them
  for(i in seq(nfolds)){
    which=foldid==i
    fitobj=sparsenet(x[!which,,drop=FALSE],y[!which],weights=weights[!which],parms=parms, ...)
    preds=predict(fitobj,x[which,,drop=FALSE])
    for(j in seq(ngamma)){
      nlami=length(fitobj$coef[[j]]$lambda)
      predmat[which,seq(nlami),j]=preds[[j]]
      nlams[i,j]=nlami
    }
    if(trace.it)cat(i)
  }
  

  N=length(y) - apply(is.na(predmat),c(2,3),sum)
  cvraw=switch(type.measure,
    "mse"=(y-predmat)^2,
    "mae"=abs(y-predmat)
    )
  
  cvm=apply(cvraw,c(2,3),weighted.mean,w=weights,na.rm=TRUE)
for(j in seq(ngamma))cvraw[,,j]=scale(cvraw[,,j],cvm[,j],FALSE)
  cvsd=apply(cvraw^2,c(2,3),weighted.mean,w=weights,na.rm=TRUE)/(N-1)
 cvsd=sqrt(cvsd)
  obj=list(lambda=t(parms[2,,]),cvm=cvm,cvsd=cvsd,cvup=cvm+cvsd,cvlo=cvm-cvsd,nzero=nz,name=type.measure,sparsenet.fit=sparsenet.object,call=this.call)
  whichmin=argmin(cvm)
  obj$parms.min=parms[,whichmin[2],whichmin[1]]
  obj$which.min=whichmin
  which=cvm< min(cvm)+cvsd
  nz[!which]=1e40
  whichcv=argmin(nz)
  obj$parms.1se=parms[,whichcv[2],whichcv[1]]
  obj$which.1se=whichcv
  class(obj)="cv.sparsenet"
  obj
}
