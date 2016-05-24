"ada.machine" <-
function(x,y,test.x,test.y,iter=50,nu=0.1,bag.frac=0.5,lossObj,oldObj=NULL,na.action=na.action,...){
  kapstat<-function(tab=diag(2) ){
    if(dim(tab)[1]==1){
      return(0)
    }
    if(dim(tab)[1]==dim(tab)[2]){
      rs<-apply(tab,2,sum)
      cs<-apply(tab,1,sum)
      N<-sum(rs)
      E<-sum(rs*cs)/N^2
      O<-sum(diag(tab))/N
      return( (O-E)/(1-E) )
    }else{
      return(0.5)
    }
  }
  tmp<-function(i){
    a1<-sample(which(y==1),1)
    a2<-sample(which(y==-1),1)
    ind<-c(a1,a2)
    return(c(sample(setdiff(1:n,ind),n-val-2,FALSE),ind))
  }

  n=dim(x)[1]
  fit=list()
  y<-as.numeric(y)
  dat<-data.frame(y=y,x)
  
  w=rep(1,n)/n
  oobm.mat<-matrix(0,nrow=n,ncol=2)
  fits=rep(0,n)
  
  atmp=alpha=vector(length=iter)
  oobm.err<-rep(0,iter)
  train.err<-rep(0,iter)
  train.kap<-rep(0,iter)
  start=0
  if(!is.null(test.y)){
    fit=oldObj$model$trees
    test.err<-rep(0,iter)
    test.kap<-rep(0,iter)
    test.n<-dim(test.x)[1]
    fits.test<-rep(0,test.n)
  }
  if(!is.null(oldObj)){
    fit=oldObj$model$trees
    w=oldObj$model$lw
    oobm.mat=oldObj$model$oob.str$oobm.mat
    fits=oldObj$model$F[[1]]
    start=oldObj$iter
    alpha[1:start]<-oldObj$model$alpha
    train.err[1:(start)]<-oldObj$model$err[,1]
    train.kap[1:(start)]<-oldObj$model$err[,2]
    oobm.err[1:(start)]<-oldObj$model$oob.str$oobm.err
    nu=oldObj$nu
    bag.frac=oldObj$bag.frac
    if(!is.null(test.y)){
      test.err[1:start]<-oldObj$model$err[,3]
      test.kap[1:start]<-oldObj$model$err[,4]
      fits.test<-oldObj$model$F[[2]]
    }
  }  
  val<-floor(bag.frac*n)
  a<-NULL
  if(val<n){
    a<-sapply(1:iter,tmp)
  } 
  if(is.vector(a)){
    a<-t(as.matrix(a))
  }
  start<-start +1
  wfun=lossObj$wfun
  coefs=lossObj$coefs
  method=lossObj$method
  predict.type=lossObj$predict.type
  shift=lossObj$shift
  f1<-f2<-0
  for (m in start:iter){
    xval=1:n
    if(!is.null(a)){
      xval=a[,m]
    }
    fit[[m]] =rpart(y~.,data=dat[xval,],weights=w[xval],method=method,x=FALSE,y=TRUE,na.action=na.action,...)
    f<-predict.type(fit[[m]],dat)
    errm=sum(w*(sign(f)!=y))
    if( (1-errm)==1 | errm==1 ){
      errm=(1-errm)*0.0001+errm*.9999
    }
    alp=0.5*log( (1-errm)/errm)
    alpha[m]=nu*coefs(w,f,y,alp)
    fits<-fits+alpha[m]*f
    if(shift){
      f1=(1-nu)*f+fits
      atmp[m]=1-nu+alpha[m]
    }
    w=wfun(y,fits)
    w=w/sum(w)
    tab<-table(sign(fits),y)
    train.err[m]<-1-sum(diag(tab))/n
    train.kap[m]<-1-kapstat(tab)
    
    indx<- setdiff(1:n,xval)
    btmp<-as.numeric(as.factor(sign(fits)[indx]))
    if(length(btmp)==1){
      oobm.mat[indx,btmp]<-oobm.mat[indx,btmp]+1
    }else{
      oobm.mat[indx,][btmp==1,1]<- oobm.mat[indx,][btmp==1,1]+1
      oobm.mat[indx,][btmp==2,2]<- oobm.mat[indx,][btmp==2,2]+1
      denom<-apply(oobm.mat,1,sum)
      vals<-denom>0
      if(sum(vals)==1){
	ytr<-c(-1,1)[which.max(oobm.mat[vals,])]
      }else{
        ytr<-c(-1,1)[apply(oobm.mat[vals,],1,which.max)]
      }
      oobm.err[m]<-sum(ytr!=y[vals])/length(vals)
    }
    if(is.null(test.y)){
      next
    }
    fit1<-predict.type(fit[[m]],test.x)
    fits.test<-fits.test +alpha[m]*fit1
    if(shift){
      f2=(1-nu)*fit1+fits.test
    }
    tab<-table(sign(fits.test),test.y)
    test.err[m]<- 1-sum(diag(tab))/test.n
    test.kap[m]<-1-kapstat(tab)
  }
  if(shift){
    alpha=atmp
    fits<-f1
    fits.test<-f2
  }
  a1=(fits==0)
  if(sum(a1)>0)
    fits[a1]<-sample(c(-1,1),sum(a1),TRUE,c(.5,.5))
  errs<-cbind(train.err,train.kap)
  ans<-list()
  ans[[1]]=fits
  if(!is.null(test.y)){
    errs<-cbind(errs,test.err,test.kap)
    ans[[2]]=fits.test
  }
  obj=list(trees=fit,alpha=alpha,F=ans,errs=errs,oob.str=list(oobm.err=oobm.err,oobm.mat=oobm.mat),lw=w,shift=shift,lossObj=lossObj)
  return(obj)
}
