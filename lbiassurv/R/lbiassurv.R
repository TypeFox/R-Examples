
lbfit.nonpar<-function(time,censor,boot=FALSE,
                       boot.control=list(quantile=TRUE,use.median=FALSE,confidence.level=0.95,iter=1000),
                       fit.control=list(iter=500,tol=1e-6)) 
{

  maxiter=fit.control$iter
  tol=fit.control$tol
  
  # check time is a numeric vector 
  # check censor is a vector of zero and ones
  # check censor has the same size as of time
  # check controls
  
  # This fitting function will be used inside bootstrap function
  fitnonpar<-function(time,censor) 
    {
    tt<-time
    d <- censor    
    #Removes zeros from a vector
    remzero<-function(x) {x[x!=0]}
    
    #Function to split time and failure indicators into two seperate vectors of
    #uncensored times and censored times, the way Vardi wrote his algorithm
    svcounts<-function(x) {
      .C("countsorted",PACKAGE="lbiassurv",as.double(x),as.integer(length(x)),y=double(length(unique(x))))$y
    }
    
    splitxy<-function(tt,delta) {
      x<-tt[delta==1]
      y<-tt[delta==0]
      list(x=x,y=y) 
    }
    
    # The algorithm starts from here
    ot<-order(tt)
    xy<-splitxy(tt[ot],d[ot])
    tj<-unique(tt[ot])
    h<-length(tj)
    m<-length(xy$x)
    n<-length(xy$y)
    iter=0
    if (n==0) {
      if (h==m) 
        vout<-list(pvec=rep(1/h,h),iter=iter)
      else
        vout<-list(pvec=svcounts(xy$x)/m,iter=iter)
    }
    else if (h==(m+n)) {
      vout<-.C("vardir",PACKAGE="lbiassurv",as.double(xy$x),as.double(xy$y),as.integer(m),as.integer(n),
               as.double(tj),as.integer(h),as.integer(maxiter), as.double(tol),
               pvec=double(h),iter=integer(1),as.double(xy$x),
               as.integer(m),as.double(xy$y),as.integer(n))
    }
    else { 
      xu<-unique(xy$x)
      yu<-unique(xy$y)
      vout<-.C("vardir",PACKAGE="lbiassurv",as.double(xy$x),as.double(xy$y),as.integer(m),as.integer(n),
               as.double(tj),as.integer(h),as.integer(maxiter), as.double(tol),
               pvec=double(h),iter=integer(1),as.double(xu),
               as.integer(length(xu)),as.double(yu),as.integer(length(yu)))
    }
    #  ret<-list(t=tj,lbpvec=vout$pvec,
    #            ubpvec=(vout$pvec/tj)/sum(vout$pvec/tj),
    #            t.break=c(0,tj),surv.break=c(1,1-cumsum((vout$pvec/tj)/sum(vout$pvec/tj))),
    #            iterations=vout$iter, conv=(iter<maxiter), method="nonpar")
#    par=cbind(tj,(vout$pvec/tj)/sum(vout$pvec/tj))
    par=cbind(tj,1-cumsum((vout$pvec/tj)/sum(vout$pvec/tj)))
    colnames(par)=c("time","surv")
    ret<-list(par=par,
              survfun=stepfun(tj,c(1,1-cumsum((vout$pvec/tj)/sum(vout$pvec/tj)))),
              iterations=vout$iter, conv=(iter<maxiter), method="nonpar")  
    class(ret)="lbsurvfit"
    return(ret)
  }
  
  bootnonpar=function(time,censor) 
    {
    B=boot.control$iter
    type=7
    conf.level=boot.control$confidence.level
    usemed=boot.control$use.median
    nonpar=boot.control$quantile
    
    a=(1-conf.level)/2
    n=length(time)
             survmat=matrix(0,nrow=B,ncol=n)
             origvardi=fitnonpar(time,censor)
             for (i in 1:B) {
               b1=sample(n,n,replace=TRUE)
               bootdata=time[b1]
               bootcens=censor[b1]
               ox1=order(bootdata)
               bootdata=bootdata[ox1]
               bootcens=bootcens[ox1]
               vardi=fitnonpar(bootdata,bootcens)
               fu1=vardi$survfun
               survmat[i,]=fu1(time)
             }
             if (nonpar) {
               lower=apply(survmat,2,quantile,probs=a,type=type)
               upper=apply(survmat,2,quantile,probs=1-a,type=type)
             }
             else {
               sdsx=apply(survmat,2,sd)
               if (!usemed) {
                 lower=pmax(0,origvardi$survfun(time)+qnorm(a)*sdsx)
                 upper=pmin(1,origvardi$survfun(time)+qnorm(1-a)*sdsx)
               }
               else {
                 medvardi=apply(survmat,2,median)
                 lower=pmax(0,medvardi+qnorm(a)*sdsx)
                 upper=pmin(1,medvardi+qnorm(1-a)*sdsx)
               }
             }
    ret=origvardi
    ot=order(time)
    ret$par=cbind(ret$par,lower[ot],upper[ot])
#    colnames(par)=c("time","surv","lower","upper")    
    ret$lowerfun=stepfun(time[ot],c(1,lower[ot]))
    ret$upperfun=stepfun(time[ot],c(1,upper[ot]))
    return(ret)
  }

if (!is.logical(boot)) stop("The boot option in bootstrap must be logical")  

  if (!(boot)) 
  {
  warning("The boot option is FALSE, no confidence bounds are produced")
  return(fitnonpar(time,censor))
  }

if (boot) 
  {
  cat("bootstrapping...")
  return(bootnonpar(time,censor))   
  }
}



lbsample <- function(n, family, par=list(shape,rate,meanlog,sdlog),censor.vec=rexp(n))
{
  shape=par$shape
  rate=par$rate
  meanlog=par$meanlog
  sdlog=par$sdlog
  censts=censor.vec
  size=n
  
  
  
  lbsample.weibull<-function(size,shape,rate,censts ) {
    g<-function (s, lambda, p) {
      (s^(1/p))/lambda
    }
    
    #length-biased time to death
    ftimes<-g(rgamma(n,shape=1+1/shape),rate,shape)
    trunctimes<-runif(n,max=ftimes)
    restimes<-ftimes-trunctimes
    #if censoring time is equal to death time, the observation is censored
    deathind<-1*(restimes<=censts)
    times<-trunctimes+deathind*restimes+(1-deathind)*censts
    ot<-order(times)
    
    #  data.frame(data=times[ot],censor=deathind[ot],onset=trunctimes[ot],study=deathind[ot]*restimes[ot]+(1-deathind[ot])*censts[ot],death=restimes[ot],censor=censts[ot])  
    # some objects are removed to be coherent with other functions
    return(list(time=times[ot],censor=deathind[ot],onset=trunctimes[ot]))  
  }
  
  
  
  lbsample.gamma<-function(size,shape,rate,censts)
  {  
    n<-size
    #might as well do this now: censoring time  
    lbtimevec=rgamma(size,shape=shape+1,rate=rate)
    lbtruncvec=runif(size,max=lbtimevec)
    residtimevec=lbtimevec-lbtruncvec
    obstimevec=lbtruncvec+pmin(residtimevec,censts)
    deltavec=1*(residtimevec<censts)
    return(list(time=obstimevec,censor=deltavec,onset=lbtruncvec))
  }
    
  lbsample.lognormal<-function(size,meanlog=0,sdlog=1,censts)
  {
    lbtimevec=exp(rnorm(size,meanlog+sdlog^2,sdlog))
    lbtruncvec=runif(size,max=lbtimevec)
    residtimevec=lbtimevec-lbtruncvec
    obstimevec=lbtruncvec+pmin(residtimevec,censts)
    deltavec=1*(residtimevec<censts)
    return(list(time=obstimevec,censor=deltavec,onset=lbtruncvec))
  }
  
  lbsample.loglogistic<-function(size,shape,rate,censts)
  { 
    alpha<-shape
    theta<-(rate)^shape
    z=rbeta(size,1-1/alpha,1+1/alpha)
    lbtimevec=((1-z)/(theta*z))^(1/alpha)
    lbtruncvec=runif(size,max=lbtimevec)
    residtimevec=lbtimevec-lbtruncvec
    obstimevec=lbtruncvec+pmin(residtimevec,censts)
    deltavec=1*(residtimevec<censts)
    return(list(time=obstimevec,censor=deltavec,onset=lbtruncvec))
  }
  

  #initial is a list of initial parameters
  if(family=="weibull"){return(lbsample.weibull(size,shape,rate,censts))}
  if(family=="gamma"){return(
    lbsample.gamma(size,shape,rate,censts))}
  
  if(family=="exponential"){return(
    lbsample.gamma(size,shape=1,rate,censts))}
  
  if(family=="lognormal"){return(
    lbsample.lognormal(size,meanlog,sdlog,censts))}
  if (family=="loglogistic"){return(
    lbsample.loglogistic(size,shape,rate,censts))}  
}


lbfit.par <- function(time, censor, family, initial=list(shape,rate,meanlog,sdlog))
{
  
  survdata=list(data=time,censor=censor)
  shape=initial$shape
  rate=initial$rate
  meanlog=initial$meanlog
  sdlog=initial$sdlog
  scale=1/initial$rate
  

  
  weibull.fit.optim=function(survdata,guess=c(0,0),...) {
    
    loglik.lb.weibull=function(log.shape=0,log.rate=0,survdata) {
      x=survdata$data
      d=survdata$censor
      llik=d*dweibull(x,shape=exp(log.shape),scale=exp(-log.rate),log=TRUE)+(1-d)*pweibull(x,shape=exp(log.shape),scale=exp(-log.rate),lower.tail=FALSE,log.p=TRUE)-lgamma(1+exp(-log.shape))+log.rate
      #attr(llik,"gradient") = c( d/shape+d*log(rate)+d*log(x)-digamma(1+1/shape)/((shape^2)*gamma(1+1/shape))-log(rate*x)*(rate*x)^shape,(1+d*(shape-1))/rate - shape*(rate^(shape-1))*x^shape)
      return(sum(llik))
    }
    
    loglik.lbweib=function(p){
      y=-loglik.lb.weibull(p[1],p[2],survdata=survdata)
      return(y) 
    }
    ret=optim(guess,loglik.lbweib,hessian=TRUE,method="BFGS")
    sdal=sqrt(diag(solve(ret$hessian)))
    retweib=list(shape.fit=exp(ret$par[1]),rate.fit=exp(ret$par[2]),shape.sd=exp(ret$par[1])*sdal[1],rate.sd=exp(ret$par[2])*sdal[2],loglike.maximum=-ret$value,hessian=ret$hessian, code=ret$convergence,iterations=ret$counts)
    return(retweib)
  }
  
  
  gamma.fit.optim=function(survdata,guess=c(0,0),...) {
    loglik.lb.gamma=function(log.shape=1,log.scale=1,survdata) {
      x=survdata$data
      d=survdata$censor
      llik=d*dgamma(x,shape=exp(log.shape),scale=exp(log.scale),log=TRUE)+(1-d)*pgamma(x,shape=exp(log.shape),scale=exp(log.scale),lower.tail=FALSE,log.p=TRUE)-log.shape-log.scale
      #attr(llik,"gradient") = c( d/shape+d*log(rate)+d*log(x)-digamma(1+1/shape)/((shape^2)*gamma(1+1/shape))-log(rate*x)*(rate*x)^shape,(1+d*(shape-1))/rate - shape*(rate^(shape-1))*x^shape)
      return(sum(llik))
    }
    
    loglik.lbgamma=function(p){
      y=-loglik.lb.gamma(p[1],p[2],survdata=survdata)
      return(y) 
    }
    #ret=nlm(loglik.lbgamma,guess,hessian=TRUE,...,survdata=survdata)
    ret=optim(guess,loglik.lbgamma,hessian=TRUE,method="BFGS")
    sdal=sqrt(diag(solve(ret$hessian)))
    retgamma=list(shape.fit=exp(ret$par[1]),scale.fit=exp(ret$par[2]),shape.sd=exp(ret$par[1])*sdal[1],scale.sd=exp(ret$par[2])*sdal[2],loglike.maximum=-ret$value,hessian=ret$hessian, code=ret$convergence,
                  iterations=ret$counts,model="gamma")
    return(retgamma)
  }
  
  
  
  exponential.fit.optim=function(survdata,guess=0,...) {
    loglik.lb.exp=function(log.rate=0,survdata) {
      x=survdata$data
      d=survdata$censor
      llik=d*dweibull(x,shape=1,scale=exp(-log.rate),log=TRUE)+(1-d)*pweibull(x,shape=1,scale=exp(-log.rate),lower.tail=FALSE,log.p=TRUE)+log.rate
      #attr(llik,"gradient") = c( d/shape+d*log(rate)+d*log(x)-digamma(1+1/shape)/((shape^2)*gamma(1+1/shape))-log(rate*x)*(rate*x)^shape,(1+d*(shape-1))/rate - shape*(rate^(shape-1))*x^shape)
      return(sum(llik))
    }
    
    loglik.lbexp=function(p){
      y=-loglik.lb.exp(p[1],survdata=survdata)
      return(y) 
    }
    #ret=nlm(loglik.lbexp,guess,hessian=TRUE,...,survdata=survdata)
    ret=optim(guess,loglik.lbexp,hessian=TRUE,method="BFGS")
    sdal=sqrt(diag(solve(ret$hessian)))
    retexp=list(rate.fit=exp(ret$par[1]),rate.sd=exp(ret$par[1])*sdal[1],loglike.maximum=-ret$value,hessian=ret$hessian, code=ret$convergence,
                iterations=ret$counts,model="exponential")
    return(retexp)
  }
  
  lognormal.fit.optim=function(survdata,guess=c(0,0),...) {
    loglik.lb.lognorm=function(mu=0,log.sigma=0,survdata) {
      x=survdata$data
      d=survdata$censor
      llik=d*dlnorm(x,meanlog=mu,sdlog=exp(log.sigma),log=TRUE)+(1-d)*plnorm(x,meanlog=mu,sdlog=exp(log.sigma),lower.tail=FALSE,log.p=TRUE)-mu-exp(2*log.sigma)/2
      #attr(llik,"gradient") = c( d/shape+d*log(rate)+d*log(x)-digamma(1+1/shape)/((shape^2)*gamma(1+1/shape))-log(rate*x)*(rate*x)^shape,(1+d*(shape-1))/rate - shape*(rate^(shape-1))*x^shape)
      return(sum(llik))
    }
    
    loglik.lblnorm=function(p){
      y=-loglik.lb.lognorm(p[1],p[2],survdata=survdata)
      return(y) 
    }
    #ret=nlm(loglik.lblnorm,guess,hessian=TRUE,...,survdata=survdata)
    ret=optim(guess,loglik.lblnorm,hessian=TRUE,method="BFGS")
    sdal=sqrt(diag(solve(ret$hessian)))
    retlnorm=list(meanlog.fit=ret$par[1],sdlog.fit=exp(ret$par[2]),meanlog.sd=sdal[1],sdlog.sd=exp(ret$par[2])*sdal[2],loglike.maximum=-ret$value,hessian=ret$hessian, code=ret$convergence,
                  iterations=ret$counts,model="lognormal")
    return(retlnorm)
  }
  
  
  
  loglogistic.fit.optim=function(survdata,guess=c(1,0),...) {
    loglik.lb.loglogis=function(log.shape=1,log.scale=0,survdata) {
      x=survdata$data
      d=survdata$censor
      llik=d*dllogis(x,shape=exp(log.shape),scale=exp(log.scale),log=TRUE)+(1-d)*pllogis(x,shape=exp(log.shape),scale=exp(log.scale),lower.tail=FALSE,log.p=TRUE)+log.shape-log.scale+log(sin(pi/exp(log.shape)))
      #attr(llik,"gradient") = c( d/shape+d*log(rate)+d*log(x)-digamma(1+1/shape)/((shape^2)*gamma(1+1/shape))-log(rate*x)*(rate*x)^shape,(1+d*(shape-1))/rate - shape*(rate^(shape-1))*x^shape)
      return(sum(llik))
    }
    
    loglik.lbllogis=function(p){
      y=-loglik.lb.loglogis(p[1],p[2],survdata=survdata)
      return(y) 
    }
    #ret=nlm(loglik.lbllogis,guess,hessian=TRUE,...,survdata=survdata)
    ret=optim(guess,loglik.lbllogis,hessian=TRUE,method="BFGS")
    sdal=sqrt(diag(solve(ret$hessian)))
    retllogis=list(shape.fit=exp(ret$par[1]),scale.fit=exp(ret$par[2]),shape.sd=exp(ret$par[1])*sdal[1],scale.sd=exp(ret$par[2])*sdal[2],loglike.maximum=-ret$value,hessian=ret$hessian, code=ret$convergence,
                   iterations=ret$counts,model="loglogistic")
    return(retllogis)
  }
  
  
  
  
  
  
  
  #initial is a list of initial parameters
  if(family=="weibull"){return(weibull.fit.optim(survdata,guess=c(log.shape=log(shape), lograte=log(rate)))) }
  if(family=="gamma"){return(gamma.fit.optim(survdata,guess=c(log.shape=log(shape), logscale=log(scale))))}
  if(family=="exponential"){return(exponential.fit.optim(survdata,guess=c(log.rate=log(rate))))}
  if(family=="lognormal"){return(lognormal.fit.optim(survdata,guess=c(mu=meanlog,log.sigma=log(sdlog))))}
  if(family=="loglogistic"){return(loglogistic.fit.optim(survdata,guess=c(shape=log(shape),scale=log(scale))))}
}


