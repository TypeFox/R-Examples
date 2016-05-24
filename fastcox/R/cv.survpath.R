######################################################################
## These functions are minor modifications or directly copied from the 
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010). 
##        Regularization Paths for Generalized Linear Models via 
##        Coordinate Descent. 
##        Journal of Statistical Software, 33(1), 1-22. 
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason they are copied here is because they are internal functions
## and hence are not exported into the global environment.
## The original comments and header are preserved.

cv.survpath=function(outlist,lambda,x,y,d,foldid){
  nfolds=max(foldid)
  cvraw=matrix(NA,nfolds,length(lambda))
  for(i in seq(nfolds)){
    which=foldid==i
    fitobj=outlist[[i]]
    coefmat=predict(fitobj,type="coeff")
    plfull=survpath.deviance(x=x,y=y,d=d,beta=coefmat)
    plminusk=survpath.deviance(x=x[!which,],y=y[!which],d=d[!which],beta=coefmat)
	cvraw[i,seq(along=plfull)]=plfull-plminusk
  }
  N=nfolds - apply(is.na(cvraw),2,sum)
  weights=as.vector(tapply(d,foldid,sum))
  cvraw=cvraw/weights
  cvm=apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
  cvsd=sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(N-1))
  list(cvm=cvm,cvsd=cvsd,name="Partial-Likelihood Deviance")
}



survpath.deviance=function(x,y,d,beta=NULL){
  storage.mode(x)="double"
  y=as.double(y)
  d=as.double(d)
  nobs=as.integer(length(y))
  nvars=as.integer(ncol(x))
  if(is.null(beta)){
    beta=double(0)
    nvec=as.integer(1)
    nvars=as.integer(0)
  }
  else{
    beta=as.matrix(beta)
    nvec=as.integer(ncol(beta))
  }
  fit=.Fortran("OBJ",nobs,nvars,x,y,d,beta,nvec,flog=double(nvec),
    PACKAGE="fastcox")
2*fit$flog*nobs
}