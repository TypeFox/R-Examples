Lifedata.MLE <-
function(formula, data, xt.dat=NULL, dist, method="BFGS", subset, truncation, na.action, weights, ref_time=NULL, starts=NULL, ...)
{

  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "truncation", "na.action", "weights"), names(mf), 0L)
  if (m[1]==0) stop("a formula argument is required")
  mf <- mf[c(1, m)]  
  mf[[1]] <- as.name("model.frame")
  mdat <- eval(mf, parent.frame())
  n=nrow(mdat)
  mt <- attr(mdat, "terms")
  
  if(is.null(xt.dat)){
    X <- model.matrix(mt, mdat)    
    ll <- ncol(X)
  } else{
    ll <- ncol(xt.dat)-1
    dist <-paste(dist, ".ce", sep="")
  }

  y <- model.extract(mdat, "response")
  truncation=model.extract(mdat, "truncation")  
  wts=model.extract(mdat, "weights")
  
  if(is.null(wts)){
    wts=rep(1, n)
  }

  if(length(wts)!=n){
    print("The length of weights is not equal to the subset")
  } 
  
  switch(dist,
         "weibull"={
           minusloglik=miniusloglik.sev.wts
           qfun=qsev
           pfun=psev
         },
         "lognormal"={
           minusloglik=miniusloglik.normal.wts
           qfun=qnorm
           pfun=pnorm
         },
         "loglogistic"={
           minusloglik=miniusloglik.logis.wts
           qfun=qlogis
           pfun=plogis
         },
         "frechet"={
           minusloglik=miniusloglik.lev.wts
           qfun=qlev 
           pfun=plev
         },
         "weibull.ce"={
           minusloglik=miniusloglik.ce.xt.sev
           qfun=qsev
         },
         "frechet.ce"={
           minusloglik=miniusloglik.ce.xt.lev
           qfun=qlev
         },
         "lognormal.ce"={
           minusloglik=miniusloglik.ce.xt.norm
           qfun=qnorm
         },
         "loglogistic.ce"={
           minusloglik=miniusloglik.ce.xt.logis
           qfun=qlogis
         })
  

  # pick start values 
  if(is.null(starts)){
    fit=survfit(y~1, weights=wts)  
    loc=kaplan.meier.location(fit)  
    coef=coef(lm(log(loc[,1])~qfun(loc[,2])))   
    starts=c(coef[1], rep(0, ll-1), log(coef[2]))
  }
  
  survdat=NULL
  if(is.null(xt.dat)){
    survdat=as.data.frame(cbind(y, wts=wts, X))
  }else{
    survdat=ce.dat.prep(xt.dat=xt.dat, failure.dat=y, ref_time=ref_time)
  }
  

  # calculate MLEs
  res.mle=lifetime.mle(dat=survdat,  minusloglik=minusloglik, starts=starts, method=method)
  
  # return values
  coef=res.mle$coef
  if(is.null(xt.dat)){
    names(coef)=c("(Intercept)", colnames(X[,-1]), "logsigma")
  } else{
    names(coef)=c("mu", colnames(xt.dat)[-c(1,2)], "logsigma")
  }
 
  rownames(res.mle$vcov)=colnames(res.mle$vcov)=NULL
  
  if(ncol(y)==2 & is.null(xt.dat)){
    zz=(log(y[,1])-as.matrix(X)%*%as.matrix(coef[1:ll]))/exp(coef[length(coef)])
    surv=as.vector(1-pfun(zz))
  }else{
    surv=NULL
  }
    
  res=list(call=mf, formula=formula, coef=coef, vcov = res.mle$vcov, min = res.mle$min, surv=surv, dat=survdat)
  class(res)="Lifedata.MLE"
  
  return(res)
}
