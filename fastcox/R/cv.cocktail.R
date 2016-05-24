cv.cocktail=function(x,y,d,lambda=NULL,nfolds=5,foldid,...){
  if(!is.null(lambda)&&length(lambda)<2)stop("Need more than one value of lambda for cv.cocktail")
  N=nrow(x)
###Fit the model once to get dimensions etc of output
  y=drop(y) # we dont like matrix responses unless we need them
  cocktail.object=cocktail(x,y,d,lambda=lambda,...)
  lambda=cocktail.object$lambda
  nz=sapply(predict(cocktail.object,type="nonzero"),length)
  if(missing(foldid)) foldid=sample(rep(seq(nfolds),length=N)) else nfolds=max(foldid)
  if(nfolds<3)stop("nfolds must be bigger than 3; nfolds=10 recommended")
   outlist=as.list(seq(nfolds))
###Now fit the nfold models and store them
  for(i in seq(nfolds)){
    which=foldid==i
    if(is.matrix(y))y_sub=y[!which,]else y_sub=y[!which]
	if(is.matrix(d))d_sub=d[!which,]else d_sub=d[!which]
    outlist[[i]]=cocktail(x[!which,,drop=FALSE],y_sub,d_sub,lambda = lambda,...)
  }
  ###What to do depends on the type.measure and the model fit
  fun=paste("cv",class(cocktail.object)[[2]],sep=".")
  cvstuff=do.call(fun,list(outlist,lambda,x,y,d,foldid))
  cvm=cvstuff$cvm
  cvsd=cvstuff$cvsd
  cvname=cvstuff$name
  
out=list(lambda=lambda,cvm=cvm,cvsd=cvsd,cvup=cvm+cvsd,cvlo=cvm-cvsd,nzero=nz,name=cvname,cocktail.fit=cocktail.object)
  lamin=getmin(lambda,cvm,cvsd)
  obj=c(out,as.list(lamin))
  class(obj)="cv.cocktail"
  obj
}
