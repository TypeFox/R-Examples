predict.TransModel<-
function(object,...){
  arg<-list(...)
  r<-object$r
  z<-object$z
  if(is.null(arg$covar)) arg$covar<-rep(0,object$nz)
  arg$covar<-as.numeric(arg$covar)
  Rt<-object$Rt

  if(is.null(arg$new.time)) arg$new.time<-object$ord.time
  n<-length(arg$new.time)
  p<-ncol(object$z)

  if(p>1){
    Rt0<-Rt*exp(-sum(object$coefficients*apply(z,2,mean)))
    #arg$covar<-arg$covar-apply(z,2,mean)
  }
  if(p==1){
    Rt0<-Rt*exp(-sum(object$coefficients*mean(z)))
    #arg$covar<-arg$covar-mean(z)
  }
  err<-Rt0*exp(sum(arg$covar*object$coefficients))
  if(r==0) St<-exp(-err)
  if(r>0) St<-exp(-1/r*log(1+r*err))
  St.fun<-stepfun(object$ord.time,c(1,St))
  new.St<-St.fun(arg$new.time)
  if(r==0) new.err<-log(-log(new.St))
  if(r>0) new.err<-log((exp(-r*log(new.St))-1)/r)

  pred<-data.frame(time=arg$new.time,survival=new.St)

  if(object$CICB.st){
    new.Ht<-matrix(ncol=n,nrow=nrow(object$rep.beta))
    raw.Ht<-matrix(ncol=length(object$ord.time),nrow=nrow(object$rep.beta))

    for(i in 1:nrow(object$rep.beta)){
      if(r>0) Rt.fun<-stepfun(object$ord.time,c(1,object$rep.Rt[i,]))
      if(r==0) Rt.fun<-stepfun(object$ord.time,c(0,object$rep.Rt[i,]))
      new.Ht[i,]<-log(Rt.fun(arg$new.time))
      raw.Ht[i,]<-log(Rt.fun(object$ord.time))
    }
    rep.bz<-matrix(object$rep.beta%*%arg$covar,ncol=n,nrow=nrow(object$rep.beta))
    new.error<-new.Ht+rep.bz
    e.sd<-apply(new.error,2,sd)
    ul.err<-new.err+qnorm(0.975)*e.sd
    ll.err<-new.err-qnorm(0.975)*e.sd
    
    raw.error<-raw.Ht+matrix(object$rep.beta%*%arg$covar,ncol=length(object$ord.time),nrow=nrow(object$rep.beta))
    Qa<-quantile(apply(abs(scale(raw.error)),1,max),0.95)
    ub.err<-new.err+Qa*e.sd
    lb.err<-new.err-Qa*e.sd
   
    if(r==0){
      ll.st<-exp(-exp(ul.err))
      ul.st<-exp(-exp(ll.err))

      lb.st<-exp(-exp(ub.err))
      ub.st<-exp(-exp(lb.err))
    } else{
      ll.st<-exp(-log(1+r*exp(ul.err))/r)
      ul.st<-exp(-log(1+r*exp(ll.err))/r)

      lb.st<-exp(-log(1+r*exp(ub.err))/r)
      ub.st<-exp(-log(1+r*exp(lb.err))/r)
    }

    ll.st[ll.st<0]<-lb.st[lb.st<0]<-0
    ul.st[ul.st>1]<-ub.st[ub.st>1]<-1
    pred<-data.frame(pred,ll.st,ul.st,lb.st,ub.st)
    pred$Qa=Qa
    }
  class(object)<-"TransModel"
  class(pred)<-"predict.TransModel"
  return(pred)
}