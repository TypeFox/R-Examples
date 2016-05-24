### lifetime.mle
## failure.threshold is the percentage
################################################################################
addt.fit=function(formula, data, initial.val=100, proc="Both", failure.threshold, time.rti=100000, method="Nelder-Mead", subset, na.action, starts=NULL,fail.thres.vec=c(70,80),...)
{
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0L)
  if (m[1]==0) stop("a formula argument is required")
  mf <- mf[c(1, m)]  
  mf[[1]] <- as.name("model.frame")
  mdat <- eval(mf, parent.frame())
  n=nrow(mdat)
  mt <- attr(mdat, "terms")
  
  #y <- model.extract(mdat, "response")  
  #X <- model.matrix(mt, mdat)    
  #ll <- ncol(X)
  
  #browser() 
  names(mdat)=c("Response", "Time", "Temp")
  dat0=as.data.frame(cbind(Temp=mdat[,"Temp"], Time=mdat[,"Time"], Response=mdat[,"Response"]))
  
  
  LS.obj=ML.obj=NULL
  
  if(proc=="LS"|proc=="Both"){
    # LSA
    LS.obj=lsa.fit(dat=dat0, initial.val=initial.val,failure.threshold=failure.threshold, time.rti=time.rti)
  }
  
  if(proc=="ML"|proc=="Both"){
    # MLA
    if(is.null(starts))
    {
      #fail.thres=50
      #ok=0
      #while(ok<1){  
      #  temp=try(lsa.fit(dat=dat0, initial.val=initial.val,failure.threshold=fail.thres, time.rti=time.rti), silent=T)
      #  ok=ifelse(class(temp)=="try-error", 0, 1)
      #  fail.thres=fail.thres+10
      #}
      #fail.thres1=fail.thres-10
      #fail.thres2=ifelse(fail.thres>100, 99, fail.thres)
      
      starts=addt.mle.initial.val(dat=dat0, time.rti=time.rti, initial.val=initial.val, failure.threshold=fail.thres.vec[1], failure.threshold1=fail.thres.vec[2]) 
    }
    ML.obj=lifetime.mle(dat=dat0, minusloglik=minus.loglik.kinetics, starts=starts, method = method , control=list(maxit=100000))
    ML.obj=mle.transform(obj=ML.obj) 
    ML.obj=c(ML.obj, rti=true.ti.compute(pars=ML.obj$coef,failure.threshold=failure.threshold,initial.val=initial.val,time.rti=time.rti)$rti)
  }
    
    
  addt.fit=list(LS.obj=LS.obj, ML.obj=ML.obj, dat=dat0, time.rti=time.rti, initial.val=initial.val,
                failure.threshold=failure.threshold)
  class(addt.fit)="addt.fit"
  
  return(addt.fit)
}

################################################################################
power.exponential.fun=function(tt, temp, alpha, beta0,beta1,gamma)
{
  x=11605/(temp+273.16)
  bb=exp(beta0+beta1*x)
  res=alpha*exp(-(tt/bb)^gamma)
  return(res)
}

################################################################################
true.ti.compute=function(pars, failure.threshold, initial.val, time.rti)
{
  pp=failure.threshold/100
  beta0=pars[2]/log(10)
  beta1=pars[3]/log(10)
  gamma=exp(pars[4])
  gamma.term=(1/gamma)*log((1-pp)/pp)/log(10)
  dK=273.16
    
  rti=beta1/(log10(time.rti)-beta0-gamma.term)-dK
  
  beta0.log10=beta0+gamma.term
  beta1.log10=beta1
  names(beta0.log10)="beta0"
  names(beta1.log10)="beta1"
  
  res=list(coef=c(beta0.log10, beta1.log10),rti=rti,time.rti=time.rti,failure.threshold=failure.threshold)
  return(res)
}

####################################################################################
print.addt.fit=function(x,...){
  obj=x
  if(!is.null(obj$LS.obj)){
    cat("Least Squares Approach: \n")
    print(round(x$LS.obj$coef0, 4))
    cat("est.TI:" , round(obj$LS.obj$rti), "\n")
  }
  
  if(!is.null(obj$ML.obj)){
    cat("\n", "Maximum Likelihood Approach:", "\n")
    
    alpha=exp(obj$ML.obj$coef[1])
    beta0=obj$ML.obj$coef[2]
    beta1=obj$ML.obj$coef[3]
    gamma=exp(obj$ML.obj$coef[4])
    
    sigma=sqrt(exp(obj$ML.obj$coef[5])^2+exp(obj$ML.obj$coef[6])^2)
    rho=exp(obj$ML.obj$coef[6])^2/(sigma^2)
    
    coef=c(alpha, beta0, beta1, gamma, sigma, rho)
    names(coef)=c("alpha", "beta0", "beta1", "gamma", "sigma", "rho")
    
    cat("Call:\n")
    print(obj$ML.obj$call)
    cat("\n")
    cat("Coefficients:\n")
    print(round(coef,4))
    cat("est.TI:", round(obj$ML.obj$rti), "\n")
    cat("\n")
    cat("Loglikelihod:\n")
    print(as.numeric(-obj$ML.obj$min))
  }
}
####################################################################################
summary.addt.fit=function(object,...)
{
  obj=object
  if(!is.null(obj$ML.obj)){
    alpha=exp(obj$ML.obj$coef[1])
    beta0=obj$ML.obj$coef[2]
    beta1=obj$ML.obj$coef[3]
    gamma=exp(obj$ML.obj$coef[4])
    
    sigma=sqrt(exp(obj$ML.obj$coef[5])^2+exp(obj$ML.obj$coef[6])^2)
    rho=exp(obj$ML.obj$coef[6])^2/(sigma^2)
    
    sigma.dev=deriv(~sqrt(exp(x)^2+exp(y)^2),c("x","y"), function(x,y){})
    rho.dev=deriv(~exp(y)^2/(exp(x)^2+exp(y)^2),c("x","y"), function(x,y){})
    
    
    ll=length(obj$ML.obj$coef)
    
    coef.dev=diag(rep(1, ll))
    coef.dev[1,1]=alpha
    coef.dev[4,4]=gamma  
    coef.dev[(ll-1),(ll-1):ll]=attr(sigma.dev(obj$ML.obj$coef[5], obj$ML.obj$coef[6]), "gradient")
    coef.dev[ll,(ll-1):ll]=attr(rho.dev(obj$ML.obj$coef[5], obj$ML.obj$coef[6]), "gradient")
    
    vcov=coef.dev%*%obj$ML.obj$vcov%*%t(coef.dev)
    
    coef=c(alpha, beta0, beta1, gamma, sigma, rho)
    names(coef)=c("alpha", "beta0", "beta1", "gamma", "sigma", "rho")
    
    mat=matrix(0, ll, 4)
    mat[,1]=coef
    mat[,2]=sqrt(diag(vcov))
    mat[,3]=mat[,1]*exp(-1.96*mat[,2]/mat[,1])
    mat[,4]=mat[,1]*exp(+1.96*mat[,2]/mat[,1])
    
    mat[c(2:3,ll),3]=mat[c(2:3,ll),1]-1.96*mat[c(2:3,ll),2]
    mat[c(2:3,ll),4]=mat[c(2:3,ll),1]+1.96*mat[c(2:3,ll),2]
    colnames(mat)=c("mean", "std", "95% Lower", "95% Upper")
    rownames(mat)=names(coef)
     
    ti.CI=addt.confint.ti.mle(obj=obj, conflevel=0.95)
    names(ti.CI)=c("est", "std", "95% Lower", "95% Upper")
    
    obj=c(obj, list(coef.mle.mat=mat), list(ti.CI=ti.CI))
  }
  
  class(obj)="summary.addt.fit"
  return(obj)    
}

####################################################################################
print.summary.addt.fit=function(x,...)
{
  if(!is.null(x$LS.obj)){
    cat("Least Squares Approach: \n")
    print(round(x$LS.obj$coef0, 4))
    cat("est.TI:" , round(x$LS.obj$rti), "\n")
    cat("Interpolation time: \n")
    print(x$LS.obj$interp.mat)
  }

  if(!is.null(x$ML.obj)){
    cat("\n", "Maximum Likelihood Approach:", "\n")
    cat("Call:\n")
    print(x$ML.obj$call)
    cat("\n")
    cat("Parameters:\n")
    print(round(x$coef.mle.mat, 4))
    cat("\n")
    cat("TI: \n")
    print(round(x$ti.CI, 4))
    cat("\n")
    cat("Loglikelihod:\n")
    print(as.numeric(-x$ML.obj$min))
  }

}

####################################################################################
plot.addt.fit=function(x, type,...){
  if(type=="data"){
    addt.data.plot(x$ML.obj$dat)
  }
  
  if(type=="ML"){
    addt.par.fitted.plot(x$ML.obj,mean.fun="kinetics",zero.correction=F,timeline.show=F,legend.pos=1)
  }
  
  if(type=="LS"){
    polynomial.interpolation(dat=x$dat,initial.val=x$initial.val,failure.threshold=x$failure.threshold)
  }
  
}

################################################################################
addt.par.fitted.plot=function(obj,mean.fun,zero.correction=F,timeline.show=T,legend.pos=1,initial.val=100,...)
{
  dat=obj$dat 
  addt.data.plot(dat=dat,zero.correction=zero.correction,timeline.show=timeline.show,legend.pos=legend.pos,...)
  
  mtime=max(dat[,"Time"])
  stime=ifelse(zero.correction,1,0)
  ww=seq(stime,mtime,,1000)
  
  coefs=obj$coef
  
  switch(mean.fun,
         "kinetics"={
           mfun=kinetics.fun
           alpha=exp(coefs[1])
           beta0=coefs[2]
           beta1=coefs[3]
           gamma=exp(coefs[4])
         },
         "kinetics0"={
           mfun=kinetics.fun
           alpha=initial.val
           beta0=coefs[1]
           beta1=coefs[2]
           gamma=exp(coefs[3])
         },
         "power.exponential"={
           mfun=power.exponential.fun
           alpha=exp(coefs[1])
           beta0=coefs[2]
           beta1=coefs[3]
           gamma=exp(coefs[4])
         },
         "power.exponential0"={
           mfun=power.exponential.fun
           alpha=initial.val
           beta0=coefs[1]
           beta1=coefs[2]
           gamma=exp(coefs[3])
         },
         
  )
  
  
  cc2=unique(dat[,"Temp"])
  
  for(i in 1:(length(cc2)))
  {
    yy=mfun(tt=ww, temp=cc2[i], alpha=alpha, beta0=beta0, beta1=beta1, gamma=gamma)
    lines(ww,yy,lwd=3,col=i)
  }
  
}
################################################################################
addt.data.plot=function(dat,zero.correction=F,timeline.show=T,legend.pos=1,xlab,ylab,...)
{
  cc=dat[,"Temp"]
  cc=as.factor(cc)
  cc1=as.numeric(cc)
  
  if(zero.correction)
  {
    dat[dat[,"Time"]==0,"Time"]=1
  }
  
  #browser()
  
  mm=max(dat[,"Response"])
  mmin=min(dat[,"Response"])
  mtime=max(dat[,"Time"])
  mintime=min(dat[,"Time"])
  
  if(!zero.correction)
  {
    mintime=0
  }
  
  plot(dat[,"Time"], dat[,"Response"],col=cc1,xlim=c(mintime,1.05*mtime),ylim=c(.8*mmin,1.2*mm),pch=cc1,xlab="Time", ylab="Responses Values",las=1,lwd=2,...)
  
  #browser()
  
  if(legend.pos==1)
  {
    legend(0.7*mtime,1.23*mm,paste(levels(cc),"C",sep=""),col=unique(cc1),pch=unique(cc1),lwd=2,lty=NA,bty="n")
  }
  
  if(legend.pos==2)
  {
    legend(mintime,.5*mm,paste(levels(cc),"C",sep=""),col=unique(cc1),pch=unique(cc1),lwd=2,lty=NA,bty="n")
  }
  
  
  if(timeline.show)
  {
    tt=unique(dat[,"Time"])
    abline(v=tt,col="grey",lty=2)
    text(tt,rep(.9*mmin,length(tt)),tt,cex=.9)
  }
}

################################################################################
lsa.fit=function(dat, initial.val, failure.threshold, time.rti)
{
  dK=273.16
  
  # polynominal fit to the data
  dat0=polynomial.interpolation(dat=dat,initial.val=initial.val,failure.threshold=failure.threshold,plot=F,plot.all=F)
  dat0=dat0[!is.na(dat0[,2]),]
  
  # least square fit to the data
  fit0=lm(I(log10(dat0[,2]))~I(1/(dat0[,1]+dK)))
  coef0=fit0$coef 
  names(coef0)=c("beta0", "beta1") 
  rti0=coef0[2]/(log10(time.rti)-coef0[1])-dK
  
  res=list(coef0=coef0, interp.mat=dat0, rti0=rti0)
  class(res)="addt.fit.lsa"
  return(res)
  
}
################################################################################
lifetime.mle=function(dat, minusloglik, starts, method = "BFGS",hessian = TRUE,...)
{
  call=match.call()
  f = function(p) {
    minusloglik(dat,p) 
  }
  oout = optim(starts, f, method = method, hessian = hessian,...)#,control=list(trace=T))
  coef = oout$par
  #browser()
  if(hessian)
  {
    vcov =solve(oout$hessian)
  }else{
    vcov=NULL
  }
  min = oout$value
  invisible(list(call = call, coef = coef,vcov = vcov, min = min,dat=dat,minusloglik=minusloglik))
}

################################################################################
# sample mean at each combo of time and temp
addt.mean.summary=function(dat)
{
  aa=tapply(dat[,3],dat[,1:2],"mean")
  
  bb=rownames(aa)
  bb=as.numeric(bb)
  cc=colnames(aa)
  cc=as.numeric(cc)
  
  mm=dim(aa)[1]
  nn=dim(aa)[2]
  
  res=data.frame(Temp=rep(bb,nn), Time=rep(cc,rep(mm,nn)), Response=as.vector(aa))
  res=res[!is.na(res[,3]),]
  res=res[order(res[,1],res[,2]),]
  rownames(res)=NULL
  return(res)
}
################################################################################
addt.confint.ti.mle=function(obj,conflevel)
{
  if(is.null(obj$ML.obj)) {
    stop("this function can only use for maximum likelihood approach")
  }
  
  failure.threshold=obj$failure.threshold
  initial.val=obj$initial.val
  time.rti=obj$time.rti
  
  obj=obj$ML.obj
  Sigma.part=obj$vcov[c(2:4),c(2:4)]
  beta0.hat=obj$coef[2]
  beta1.hat=obj$coef[3]
  pp=failure.threshold/100
  gamma.hat=exp(obj$coef[4])
  gamma.term=(1/gamma.hat)*log((1-pp)/pp)
  ti.hat=true.ti.compute(pars=obj$coef,failure.threshold=failure.threshold,initial.val=initial.val,time.rti=time.rti)$rti
  
  beta0.log10=beta0.hat/log(10)
  beta1.log10=beta1.hat/log(10)
  gamma.log10=gamma.hat*log(10)
  gamma.term.log10=log((1-pp)/pp)/gamma.log10
  
  den=log10(time.rti)-beta0.log10-gamma.term.log10
  partial.b0=beta1.log10/(den^2)
  partial.b1=1/den
  partial.gamma=-(beta1.log10/(den^2))*(gamma.term.log10/gamma.log10)
  
  Sigma.part.trans=diag(c(1/log(10),1/log(10),log(10)*gamma.hat))%*%Sigma.part%*%diag(c(1/log(10),1/log(10),log(10)*gamma.hat))
    
  Sigma.ti=t(c(partial.b0,partial.b1,partial.gamma))%*%Sigma.part.trans%*%c(partial.b0,partial.b1,partial.gamma)
  CI=c(ti.hat-qnorm(0.5+conflevel/2)*sqrt(Sigma.ti),ti.hat+qnorm(0.5+conflevel/2)*sqrt(Sigma.ti))
  #cat(100*(1-conflevel), "% CI is ", "(", round(CI[1], 3), ",", round(CI[2], 3), ") \n", sep="")
  
  CI[1]=ifelse(CI[1]<0, 0, CI[1])
  
  res=c(ti.hat, sqrt(Sigma.ti), CI)
  names(res)=c("est.", "s.e.", "lower", "upper")
  return(res)
}

################################################################################
addt.predint.ybar.mle=function(obj,conflevel,num.fut.obs=5,temp,tt)
{
  if(is.null(obj$ML.obj)) {
    stop("this function can only use for maximum likelihood approach")
  }
  
  obj=obj$ML.obj
  alpha.hat=exp(obj$coef[1])
  beta0.hat=obj$coef[2]
  beta1.hat=obj$coef[3]
  gamma.hat=exp(obj$coef[4])
  
  sigma.sq.hat=exp(obj$coef[5])^2+exp(obj$coef[6])^2
  rho.hat=exp(obj$coef[6])^2/sigma.sq.hat
  
  Sigma.betaf=diag(c(alpha.hat,1,1,gamma.hat))%*%obj$vcov[c(1:4),c(1:4)]%*%diag(c(alpha.hat,1,1,gamma.hat))
  
  ybar.hat=kinetics.fun(tt=tt, temp=temp, alpha=alpha.hat, beta0=beta0.hat, beta1=beta1.hat, gamma=gamma.hat)
  partial.vec=first.dev.kinetics.fun(tt=tt, temp=temp, alpha=alpha.hat, beta0=beta0.hat, beta1=beta1.hat, gamma=gamma.hat)
  
  Sigma.ybar=sigma.sq.hat*(rho.hat+(1-rho.hat)/num.fut.obs)+
    t(partial.vec)%*%Sigma.betaf%*%partial.vec
  
  return(c(ybar.hat-qnorm(0.5+conflevel/2)*sqrt(Sigma.ybar),ybar.hat+qnorm(0.5+conflevel/2)*sqrt(Sigma.ybar)))  
}


################################################################################
first.dev.kinetics.fun=function(tt, temp, alpha, beta0, beta1, gamma){
  x=1/(temp+273.16)
  
  partial.beta0 = alpha*gamma*(tt^gamma)*
    (1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-2)*exp(-(beta0+beta1*x)*gamma)
  partial.beta1 = alpha*x*gamma*(tt^gamma)*
    (1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-2)*exp(-(beta0+beta1*x)*gamma)
  partial.gamma = -alpha*(tt^gamma)*
    (1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-2)*exp(-(beta0+beta1*x)*gamma)*
    (log(tt)-(beta0+beta1*x))
  partial.alpha=(1+(tt^gamma)*exp(-(beta0+beta1*x)*gamma))^(-1)
  
  return(c(partial.alpha, partial.beta0, partial.beta1, partial.gamma))
}


################################################################################
# initial.val is specified if no records available at time 0
addt.data.normalization=function(dat, initial.val)
{
  temps=unique(dat[,"Temp"])
  nn=length(temps)
  
  if(!any(dat[,"Time"]==0))
  {
    tmp1=cbind(Temp=temps, Time=0, Response=100)
    dat[,"Response"]=dat[,"Response"]/initial.val*100
    dat=rbind(tmp1,dat)
    dat=dat[order(dat[,"Temp"],dat[,"Time"]),]
    rownames(dat)=NULL
    res=dat
  }
  
  if(any(dat[,"Time"]==0))
  {
    initial.val=dat[dat[,"Temp"]==min(temps) & dat[,"Time"]==0,"Response"]
    
    res=NULL
    for(i in 1:nn)
    {
      xtmp=dat[dat[,"Temp"]==temps[i],]
      
      if(any(xtmp[,"Time"]==0))
      {
        xtmp[,"Response"]=xtmp[,"Response"]/xtmp[xtmp[,"Time"]==0,"Response"]*100
      }else{
        xtmp[,"Response"]=xtmp[,"Response"]/initial.val*100
        xtmp=rbind(xtmp,c(temps[i],0,100))
      }
      res=rbind(res,xtmp)
    }
    
    res=res[order(res[,"Temp"],res[,"Time"]),]
    rownames(res)=NULL
  }
  return(res)
}


################################################################################
polynomial.interpolation=function(dat,initial.val=100,failure.threshold=80,plot.all=T,plot=T)
{
  dat=addt.mean.summary(dat=dat)
  dat=addt.data.normalization(dat=dat, initial.val=initial.val)
  
  temps=unique(dat[,"Temp"])
  nn=length(temps)
  tres=rep(NA,nn)
  
  if(plot)
  {
    if(plot.all)
    {
      plot(0,0,type="n",xlim=c(0,max(dat[,"Time"])),ylim=c(0,100),ylab="Response (Relative %)",xlab="Time (hours)")
      abline(h=failure.threshold)
    }else{
      par(mfrow=c(ceiling(nn/round(sqrt(nn))),round(sqrt(nn))))
    }
  }
  
  for(i in 1:nn)
  {
    idx=(dat[,"Temp"]==temps[i])
    
    if(sum(idx)>=3)
    {
      yy=dat[idx,"Response"]
      xx=dat[idx,"Time"]
      #browser()
      if(sum(idx)==3 & min(yy)<failure.threshold)
      {
        afit=lm(yy~I(xx)+I(xx^2))
        coefs=afit$coef
        
        tt=seq(0,max(xx),,1000)
        yyhat=coefs[1]+coefs[2]*tt+coefs[3]*tt^2
        
        
        ctmp=polyroot(c(coefs[1]-failure.threshold,coefs[2:3]))
        
        #print(ctmp)
        ctmp1=Re(ctmp[abs(Im(ctmp))<1e-5])
        ctmp2=ctmp1[ctmp1>0 & ctmp1<=max(xx)]
        if(length(ctmp2)>0)
        {
          tres[i]=min(ctmp2)
        }else{
          ff=approx(yyhat, tt, failure.threshold)
          tres[i]=min(ff$y)
        }
        #browser()
        
        if(plot)
        {
          if(!plot.all)
          {
            plot(tt,yyhat,type="l",ylim=c(0,100),xlim=c(0,max(xx)),xlab="Time",ylab="Response (Relative %)")
            abline(h=failure.threshold)
            points(xx,yy)
          }else{
            lines(tt,yyhat,lwd=2,col=i)
            points(xx,yy,col=i,pch=i,lwd=2)
          }
          points(tres[i],failure.threshold,pch=13,col=2,lwd=2)
        }
      }
      
      if(sum(idx)>3 & min(yy)<failure.threshold)
      {
        afit=lm(yy~I(xx)+I(xx^2)+I(xx^3))
        coefs=afit$coef
        
        tt=seq(0,max(xx),,1000)
        yyhat=coefs[1]+coefs[2]*tt+coefs[3]*tt^2+coefs[4]*tt^3
        #ff=approx(yyhat, tt, failure.threshold)
        #tres[i]=min(ff$y)
        #browser()
        ctmp=polyroot(c(coefs[1]-failure.threshold,coefs[2:4]))
        #print(ctmp)
        ctmp1=Re(ctmp[abs(Im(ctmp))<1e-5])
        ctmp2=ctmp1[ctmp1>0 & ctmp1<=max(xx)]
        if(length(ctmp2)>0)
        {
          tres[i]=min(ctmp2)
        }else{
          ff=approx(yyhat, tt, failure.threshold)
          tres[i]=min(ff$y)
        }
        
        
        if(plot)
        {
          if(!plot.all)
          {
            plot(tt,yyhat,type="l",ylim=c(0,100),xlim=c(0,max(xx)),xlab="Time",ylab="Response (Relative %)")
            abline(h=failure.threshold)
            points(xx,yy)
          }else{
            lines(tt,yyhat,lwd=2,col=i)
            points(xx,yy,col=i,pch=i,lwd=2)
          }
          points(tres[i],failure.threshold,pch=13,col=2,lwd=2)
        }
        
      }
    }
  }
  
  res=cbind(temps,tres)
  colnames(res)=c("Temp","Time")
  
  if(plot.all)
  {
    legend(0.7*max(dat[,"Time"]),100, paste("T_",res[,"Temp"],"=",round(res[,"Time"]),sep="") ,pch=1:nn,col=1:nn,lty=1,lwd=2,bty="n")
  }
  return(res)
  
}

################################################################################
addt.mle.initial.val=function(dat=dat, time.rti=100000, initial.val=100, failure.threshold=70, failure.threshold1=50)
{
  #browser()
  base.factor=log10(exp(1))
  dK=273.16
  
  yy=dat[,"Response"]
  tt=dat[,"Time"]
  temp=dat[,"Temp"]
  
  #initial value for alpha
  alpha=mean(yy[tt==0])
  #initial value for beta0, beta1 and gamma
  tmp.sfit1=lsa.fit(dat=dat,time.rti=time.rti,initial.val=initial.val,failure.threshold=failure.threshold)
  tmp.sfit2=lsa.fit(dat=dat,time.rti=time.rti,initial.val=initial.val,failure.threshold=failure.threshold1)
  
  pp1=failure.threshold/100
  pp2=failure.threshold1/100
  
  coef1=tmp.sfit1$coef0
  coef2=tmp.sfit2$coef0
  
  #initial value for beta1
  beta1=coef1[2]
  sbeta0=log10(time.rti)-beta1/(tmp.sfit2$rti0+dK)
  
  cc1=log((1-pp1)/pp1)
  cc2=log((1-pp2)/pp2)
  kk1=coef1[1]/base.factor
  kk2=sbeta0/base.factor
  
  kk3=(kk1-kk2)/(cc1-cc2)
  gamma=1/kk3
  
  if(gamma<0)
  {
    gamma=1
  }
  
  beta0=kk1-kk3*cc1
  beta1=beta1/11605/base.factor
  beta0=beta0+beta1*11605/(min(temp)+dK)
  
  aa=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  des=yy-aa
  sigma=sqrt(mean(des^2))
  
  start.val0=c(log(alpha), beta0, beta1, log(gamma), log(sigma))
  names(start.val0)=NULL
  
  #browser()
  
  tmp.fit=lifetime.mle(dat=dat, minusloglik=minus.loglik.kinetics.no.cor, starts=start.val0, method = "Nelder-Mead", control=list(maxit=10000))
  
  
  #cat("-log likelihood without correlation:", tmp.fit$min, "\n")
  
  tcoef=tmp.fit$coef
  alpha=exp(tcoef[1])
  beta0=tcoef[2]
  beta1=tcoef[3]
  gamma=exp(tcoef[4])
  
  aa1=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  des1=yy-aa1
  xtmp=tapply(des1,paste(temp,tt),"mean")
  sigma1=sd(xtmp)
  
  res=c(tcoef,log(sigma1))
  
  return(res)
}
################################################################################
minus.loglik.kinetics.no.cor=function(dat,pars)
{
  #print(pars)
  yy=dat[,"Response"]
  tt=dat[,"Time"]
  temp=dat[,"Temp"]
  
  alpha=exp(pars[1])
  beta0=pars[2]
  beta1=pars[3]
  gamma=exp(pars[4])
  sigma=exp(pars[5])
  sigma1=0#exp(pars[6]) #batch level
  ####
  
  #yy=yy/100
  #A=A/100
  
  #######
  aa=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  time.rti=yy-aa
  
  #plot(tt,aa)
  
  #browser()
  
  temp.time=paste(dat[,"Temp"],"C",dat[,"Time"],"H",sep="")
  cc=unique(temp.time)
  nn=length(cc)
  
  ###
  res=0
  for(i in 1:nn)
  {
    #browser()
    idx=(temp.time==cc[i])
    pp=sum(idx)
    SS=diag(sigma^2,pp)+sigma1^2*matrix(1,pp,pp)*as.numeric(pp>1)
    SS.inv=solve(SS)
    tdd=time.rti[idx]
    tres=-0.5*pp*log(2*pi)-0.5*log(det(SS))-0.5*t(tdd)%*%SS.inv%*%tdd
    #print(tres)
    res=res+tres
    #print(round(SS.inv,3))
    #browser()
  }
  
  res=as.vector(res)
  res=(-1)*res
  #print(res)
  return(res)
}


################################################################################
minus.loglik.kinetics=function(dat,pars)
{
  #print(pars)
  yy=dat[,"Response"]
  tt=dat[,"Time"]
  temp=dat[,"Temp"]
  
  alpha=exp(pars[1])
  beta0=pars[2]
  beta1=pars[3]
  gamma=exp(pars[4])
  sigma=exp(pars[5])
  sigma1=exp(pars[6]) #batch level
  ####
  
  #yy=yy/100
  #A=A/100
  
  #######
  aa=kinetics.fun1(tt=tt, temp=temp, alpha=alpha,  beta0=beta0, beta1=beta1, gamma=gamma)
  time.rti=yy-aa
  
  #plot(tt,aa)
  
  #browser()
  
  temp.time=paste(dat[,"Temp"],"C",dat[,"Time"],"H",sep="")
  cc=unique(temp.time)
  nn=length(cc)
  
  ###
  res=0
  for(i in 1:nn)
  {
    #browser()
    idx=(temp.time==cc[i])
    pp=sum(idx)
    SS=diag(sigma^2,pp)+sigma1^2*matrix(1,pp,pp)*as.numeric(pp>1)
    SS.inv=solve(SS)
    tdd=time.rti[idx]
    tres=-0.5*pp*log(2*pi)-0.5*log(det(SS))-0.5*t(tdd)%*%SS.inv%*%tdd
    #print(tres)
    res=res+tres
    #print(round(SS.inv,3))
    #browser()
  }
  
  res=as.vector(res)
  res=(-1)*res
  #print(res)
  return(res)
}

################################################################################
mle.transform=function(obj)
{
  dat=obj$dat
  coefs=obj$coef
  vcovs=obj$vcov
  
  temp=dat[,"Temp"]
  x=11605/(temp+273.16)
  mx=min(x)
  
  mat=diag(length(coefs))
  mat[2,3]=-mx
  mat[3,3]=11605
  coefs=mat%*%coefs
  coefs=as.vector(coefs)
  vcovs=mat%*%vcovs%*%t(mat)
  
  coefs->obj$coef
  vcovs->obj$vcov
  
  return(obj)
}

################################################################################
kinetics.fun=function(tt, temp, alpha, beta0, beta1, gamma)
{
  
  #browser()
  x=1/(temp+273.16)
  
  mu=beta0+beta1*x
  nu=1/gamma
  zz=(log(tt)-mu)/nu
  res=alpha/(1+exp(zz))
  return(res)
}

################################################################################
#used for MLE
kinetics.fun1=function(tt, temp, alpha, beta0, beta1, gamma)
{
  
  #browser()
  x=11605/(temp+273.16)
  mx=min(x)
  mu=beta0+beta1*(x-mx)
  nu=1/gamma
  zz=(log(tt)-mu)/nu
  res=alpha/(1+exp(zz))
  return(res)
}