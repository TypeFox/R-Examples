  # bspline smoothing function
  # x = predictor
  # degree = 1 for polygons
  # nk = knots in the range of x including the min(x) and max(x)
  require(splines)
  require(survival); require(Epi)
  
  bs.design= function(x,nk,degree){
    xl= min(x)-0.00001
    xr = max(x)+0.00001
    dx = (xr-xl)/(nk-1)
    knots = seq(xl- degree*dx,xr+degree*dx,by=dx )
    bs = spline.des(knots,x,degree+1)$design
    return(bs)
  }
  
  
  # .................................................................
  # ...........................................    Poisson smoothing:  (with possible covariates)
  # .................................................................
  #
  # v = starting value
  bspois.basic = function(x,y,offset=0,v= NULL, degree=1,nk=31,lambda=10, phi=1)
  {
  #if with covariates
  x<-as.matrix(x)
  ncov<-ncol(x)-1
  
  # Design matrix
  Z = bs.design(x[,1],nk,degree)
  
  X = as.matrix(x[,-1]) #linear covariate  (centered)    
  
  
  # second difference
  D = diag(ncol(Z))
  for (i in 1:2) D = diff(D)
  Rinv = t(D) %*% D
  
  # starting
  m0 = mean(y*exp(-offset))
  #add fixed effect
  beta = as.matrix(rep(c(0.0000),ncov))
    
  if (is.null(v)) v = rep(log(m0), ncol(Z))
  
  for (i in 1:25){
    hmu = offset+ c(Z %*% v)+ c(X %*% beta) #added F.E.
    mu = exp(hmu)
    W = mu
    Y = hmu + (y-mu)/mu
  
    # update R.E.
    ZWZ = t(Z*W)%*% Z
    ZWZRi = solve(ZWZ+ lambda*Rinv)
    v =  ZWZRi %*% t(Z*W)%*%(Y-offset - X %*% beta)
          #inv(Z'WZ+lambda*Rinv)((ZW)'(Y-offset-Xbeta)
    # update F.E.
    if (ncov>0) beta = solve(t(X*W) %*% X) %*% (t(X*W) %*% (Y - offset- Z %*% v))
          #inv(X'WX)(X'W(Y-Zv))  
   }
    hmu = offset + c(Z %*% v)   
    mu = exp(hmu)
  df = sum(ZWZRi * ZWZ) # = tr(ZWZRi %*% t(ZWZ)) = tr(ZWZRi %*% ZWZ)
  tr = -sum(log(eigen(ZWZRi)$val))
  dH= apply((Z %*% ZWZRi)*Z, 1, sum)  ## diag(Z %*% ZWZRi %*% t(Z))
  if (ncov==0) beta = NULL
  return(list(fit=mu, linear.predictor=hmu,beta = beta, v=v, df=df, dH=dH, tr=tr,phi=phi))
  }
  
  #.................................................................................
  # ............................automatic Poisson smoothing, allowing overdispersion (con covariate)
  #.................................................................................
  
  bspois= function(x,y,offset=0,nk=31,degree=1,lambda=NULL,phi=NULL,
                   alpha=0.05, err=0.0001,verbose=T,l0=l0)
  {
  
  # if lambda is provided - but phi is not:  phi is assumed 1
  if(!is.null(lambda)){ 
    if(is.null(phi)) phi=1
    fit= bspois.basic(x,y,offset,nk=nk,degree=degree,lambda=lambda,phi=phi)
    dev = 2*(log(dpois(y,y))- log(dpois(y,fit$fit)))
    # if (verbose){cat('phi=',phi)} 

  }
  
  
  # .................................  if lambda is NULL, but phi is given
  if(is.null(lambda) & !is.null(phi))
  {  
   nx = length(x)
   if (verbose){cat('Iterations: relative error in sv2-hat =',err,'\n')}
  
  # start
  lambda=l0; oldval = 2*lambda
  while (abs(oldval-lambda)/lambda > err & lambda< 1/err){
   oldval = lambda
   fit = bspois.basic(x,y,offset,nk=nk,degree=degree, lambda=lambda, phi=phi)
   df = fit$df
   tr = fit$tr 
   v = fit$v; vRv = sum(diff(diff(v))^2)
   dev = 2*(log(dpois(y,y))- log(dpois(y,fit$fit)))
   sv2 = vRv/(df-2)
    lambda= phi/sv2; lambda   
 
  if (verbose){cat('sv2=',sv2,'  df=',df,'  lambda=',lambda,'\n')}
  } # end while
  } # end if lambda=NULL
  
  
  # ........................   if both lambda and phi are to be computed
  if(is.null(lambda) & is.null(phi)){   
  nx = length(x)
  if (verbose){cat('Iterations: relative error in phi-hat =',err,'\n')}
  
  # start
  phi=1; oldphi=2*phi
   lambda=10 ; oldval = 2*lambda
 while ( (abs(oldphi-phi)/phi > err | abs(oldval-lambda)/lambda > err) & lambda< 1/err ){  #
   oldval = lambda
   oldphi = phi 
   fit = bspois.basic(x,y, offset,nk=nk,degree=degree, lambda=lambda, phi=phi)
   df = fit$df
   tr = fit$tr 
   v = fit$v; vRv = sum(diff(diff(v))^2)
   dev = 2*(log(dpois(y,y))- log(dpois(y,fit$fit)))
    phi = var( (y-fit$fit)/sqrt(fit$fit), na.rm=T)
   
   sv2 = vRv/(df-2)
    lambda= phi/sv2; lambda       #phi/sv2; 
    if (verbose){cat('phi=',phi, '  sv2=',sv2,'  df=',df,'  lambda=',lambda,'\n')}
   } # end while
  } # end if lambda=NULL
  
  # 100(1-alpha) CI for the fit
  up = fit$lin + qnorm(1-alpha/2)*sqrt(fit$dH); up = exp(up)
  lo = fit$lin - qnorm(1-alpha/2)*sqrt(fit$dH); lo = exp(lo)
  sv2=1/lambda  # solo nel 1 caso ?
  
  return(list(fit=fit$fit, linear.pred = fit$lin, 
         beta=fit$beta,v=fit$v, phi=phi, sv2=sv2, df=fit$df, 
         dH=fit$dH, lower.ci = lo, upper.ci = up, deviance= sum(dev)))
  
  } # end function
  
  
  bshazard_fit<-function(mf,covs,nbin=NULL, nk=31,degree=1,l0=10,lambda=NULL,phi=NULL,
                   alpha=0.05, err=0.0001,verbose=T)       
  {  # UseMethod("bshz")
      #   mf<-model.frame(formula=formula,data=data)
      #    covs<-  model.matrix(attr(mf,"terms"),data=mf)
          names.cov<-make.names(dimnames(covs)[2][[1]][-1])
          ncov<-(ncol(covs)-1) #number of covariate in model
       #checks of input
       if (attr(mf$Surv, "type") %in% c("interval","interval2")) print("Interval censored data not considered")
         nobs<- nrow(mf$Surv)
       #  save the Surv object and prepare it to use it in Lexis macro
            surv<- mf$Surv   
            if (ncol(surv)==2) {entry_<-rep(0,length(surv[,"time"])); exit_<-surv[,"time"]}
            if (ncol(surv)==3) {entry_<-surv[,"start"]; exit_<-surv[,"stop"]       }
            dataf<-data.frame(entry_=entry_,exit_=exit_,status=surv[,"status"],covs)
      #  compute the means of each covariate and center the covariates for estimation
            medie<-colMeans(dataf)[names.cov]  
          if (length(names.cov)>0)  dataf[,paste(names.cov,"_c",sep="")]<- dataf[,names.cov] -  matrix(medie,nrow(covs),ncol(covs)-1,byrow=T)
     
     ############# if not covariates & not pre bin ###########################
      if (length(names.cov)==0 & is.null(nbin))  {
                      fit = survfit(Surv(dataf$entry_,dataf$exit_,dataf$status)~1)
                      dt<-diff(c(0,fit$time))
                      X = fit$time
                      y = fit$n.event
                      PY = fit$n.risk* dt   # person-years at risk
                      rate = y/PY
       input<-data.frame(X=X,n.event=y,person.time=PY,raw.hazard=rate)
       input<-input[input$X>0 & input$person.time>0,]
       X<-input$X ;       y<-input$n.event    ;    PY<-input$person.time
            }     else { ###    if covariates
     ############# with covariates or prebin ###########################
     
               #  Lexis object                                                          
                  lex<-Lexis(entry=list(per=entry_),exit=list(per=exit_),exit.status = status,data=dataf)
               #  breaks to split the data according to min,max and number of bins
              
                   if (!is.null(nbin)) xx<- seq(min(entry_),max( exit_),len=(nbin+1)) else {
                    # if nbin is not specified split the data at each time at exit in the data
                   xx<- sort(unique(dataf$exit)) # [dataf$status==1] 
                    }
               #  split according to bins/times  
                   splitted<-splitLexis(lex, breaks=c(min(entry_),xx,max( exit_)),time.scale="per") 
                   splitted$per.mid<-timeBand(splitted, time.scale="per", type="middle") #take the middle of each bins as ref time
                
                # Aggregate data 
                if (length(names.cov)>0) {  # if account for covariates
                  ysum_<-aggregate(x=as.numeric(status(splitted,at="exit")==1),by=splitted[,c("per.mid",paste(names.cov,"_c",sep=""))], FUN=sum)
                  PY_  <-aggregate(x=splitted$lex.dur,                         by=splitted[,c("per.mid",paste(names.cov))], FUN=sum) 
                  } else   {               #  if no covariates 
                  ysum_<-aggregate(x=as.numeric(status(splitted,at="exit")==1),by=list(splitted[,"per.mid"]), FUN=sum)
                  PY_  <-aggregate(x=splitted$lex.dur,                         by=list(splitted[,"per.mid"]), FUN=sum)
                  }       
                  
                  ysum<-ysum_[,ncol(ysum_)]  #last column is the sum of events in each bin
                  PY<-PY_[,ncol(PY_)]        #last column is the sum of person time in each bin
                  X<-as.matrix(ysum_[,-ncol(ysum_)])   #time and  CENTERED covariates values  (if present)
                 # Order data by time 
                   ord= order(X[,1])
                   X<-X[ord,];y<-ysum[ord];PY<-PY[ord]
                   rate<-y/PY
                  input = data.frame(X,n.event=y,person.time=PY,raw.hazard=rate)
                }
       fit =  bspois(X,y,offset = log(PY),  nk=nk, degree=degree,lambda=lambda,phi=phi,
                         alpha=alpha, err=err,verbose=verbose,l0=l0)
    #output 
      #1.raw.data      
          raw.data<-input
          names(raw.data)[1] <- "time"
         
          if (ncol(covs)>1) {   #if there are covariates attach the value of covaiates
          raw.data[,paste(names.cov)]<-raw.data[,paste(names.cov,"_c",sep="")]+matrix(medie,nrow=nrow(raw.data),ncol=ncol(covs)-1,byrow=T)
          col.d<-seq(2,2+ncol(covs)-2,1)
          raw.data<-raw.data[,-col.d]
          }
      #2. fitted smoothing and parameters
        hazard=fit$fit/PY; lower.ci=fit$low/PY;   upper.ci=fit$upper/PY 
        dH<-fit$dH
         times<-as.matrix(X)[,1]
   
      if (length(names.cov)>0)  {
          cc<-medie
              #d.unique<-unique(data.frame(times,hazard,lower.ci,upper.ci,dH))
              times.unique<-times[!duplicated(times)];
              hazard.unique<-hazard[!duplicated(times)]
              lower.ci.unique<-lower.ci[!duplicated(times)]
              upper.ci.unique<-upper.ci[!duplicated(times)]
              dH.unique<-dH[!duplicated(times)]
               }   else       {
                  times.unique<-times  ;hazard.unique<-hazard;  lower.ci.unique<-lower.ci  ;upper.ci.unique<-upper.ci ;dH.unique<-dH
                  cc<-NULL}
#   return(list(raw.data=raw.data,
#   time=times.unique,cov.value=cc,hazard=hazard.unique, lower.ci=lower.ci.unique, upper.ci=upper.ci.unique
#                ,dH=dH.unique,fit=fit$fit,  beta=fit$beta, sv2=fit$sv2,phi=fit$phi, df=fit$df
       
   return(list(raw.data=raw.data,  n=nobs,
   time=times.unique,cov.value=cc,hazard=hazard.unique, lower.ci=lower.ci.unique, upper.ci=upper.ci.unique
                ,dH=dH.unique,fitted.values=fit$fit,coefficients=fit$beta, sv2=fit$sv2,phi=fit$phi, df=fit$df
           )) 
  }
############ classes and methods #####################
bshazard<-function(formula,data = parent.frame(),nbin=NULL, nk=31,degree=1,l0=10,lambda=NULL,phi=NULL,
                             alpha=0.05, err=0.0001,verbose=T) UseMethod("bshazard")
 

bshazard.default<-function(formula,data = parent.frame(),nbin=NULL, nk=31,degree=1,l0=10,lambda=NULL,phi=NULL,
                             alpha=0.05, err=0.0001,verbose=T) {
       mf<-model.frame(formula=formula,data=data)
        covs<-  model.matrix(attr(mf,"terms"),data=mf)
        y<-model.response(mf)
  est<-bshazard_fit(mf,covs,nbin, nk,degree=degree,l0,lambda,phi,alpha, err,verbose)       
  est$call<-match.call()
  class(est)<-"bshazard"
  est
}
 


plot.bshazard<-function(x,conf.int=T,overall=T,col=1,lwd=1,xlab="Time", ylab="Hazard rate",...){
 if (overall==T){
  plot(x$time,x$hazard,xlab=xlab,type="l", ylab=ylab,lwd=lwd,lty=1, col=col)
    if (conf.int==T) {
        lines(x$time, x$low, lty=2, col=col,lwd=lwd)
        lines(x$time, x$up, lty=2, col=col,lwd=lwd)
    }#se CI
 }#se overall
 if (overall==F & !is.null(x$cov.value)) {
 covs<-unique(x$raw.data[,attr(x$cov.value,"names")])
    lin<-(covs-matrix(as.numeric(x$cov.value),nrow(covs),ncol(covs),byrow=T))*matrix(x$coefficients,nrow(covs),ncol(covs),byrow=T)
    HR<-exp(rowSums(lin))
  plot(x$time,x$hazard,xlab=xlab,type="n", ylab=ylab,lwd=lwd,lty=1, col=col)
    for (i in 1:nrow(covs)) {
          h<-x$hazard*HR[i]
          lines(x$time,h,xlab=xlab,type="l", ylab=ylab,lwd=lwd,lty=1, col=col)
  }
 }
}

lines.bshazard<-function(x,conf.int=T,overall=T,col=1,lwd=1,...){
 if (overall==T){
  lines(x$time,x$hazard,lwd=lwd,lty=1, col=col)
    if (conf.int==T) {
        lines(x$time, x$low, lty=2, col=col,lwd=lwd)
        lines(x$time, x$up, lty=2, col=col,lwd=lwd)
    }#se CI
 }#se overall
 if (overall==F & !is.null(x$cov.value)) {
 covs<-unique(x$raw.data[,attr(x$cov.value,"names")])
    lin<-(covs-matrix(as.numeric(x$cov.value),nrow(covs),ncol(covs),byrow=T))*matrix(x$coefficients,nrow(covs),ncol(covs),byrow=T)
    HR<-exp(rowSums(lin))
# lines(x$time,x$hazard,xlab=xlab,type="n", ylab=ylab,lwd=lwd,lty=1, col=col)
    for (i in 1:nrow(covs)) {
          h<-x$hazard*HR[i]
          lines(x$time,h,lwd=lwd,lty=1, col=col)
  }
 }
}

print.bshazard<-function(x,...){
#  cat("Call:\n")
   cat("Call: ")
   print(x$call)
#  dput(x)
  cat("\n")
  cat("Number of records: ",x$n)
   cat("\n")
   ncov<-length(x$cov.value)
  if (is.null(x$cov.value)) {
      #n<-length(x$raw.data$time)
     n.e<-sum(x$raw.data$n.event)
     pt<-sum(x$raw.data$person.time)
     r<-n.e/pt
     print(data.frame("n.events"=as.numeric(n.e),"person.time"=as.numeric(pt),"rate"=as.numeric(r)))
  }
  if (!is.null(x$cov.value)){
    # n.e<-tapply(x$raw.data$n.event,x$raw.data[,5:(5+ncov-1)],sum)
    # pt<-tapply(x$raw.data$person.time,x$raw.data[,5],sum)
     agg<-aggregate(x$raw.data, by = list(x$raw.data[,5:(5+ncov-1)]), FUN = "sum")
     agg$rate<-agg$n.event/agg$person.time
     agg[,attr(x$cov.value,"names")]<-agg[,1:(1+ncov-1)]
     print(data.frame(agg[,c(attr(x$cov.value,"names"),"n.event","person.time","rate")]))
  } 
   # cat("Number of records",as.numeric(n),sep=" ")
#    cat("\n")
#   cat("Number of events", as.numeric(n.e),sep=" ");cat("\n")
#   cat("Total person-time", as.numeric(pt),sep=" ");cat("\n")
#    cat("Overall rate of event",as.numeric(r),sep=" ") ;cat("\n")                                    
 }

summary.bshazard<-function(object,digits=4,...){
  #  cat("Call:\n")
#  cat("Call: ")
#  print(x$call)
  #  dput(x)
        #x$raw.data$timex$raw.data$n.event,x$raw.data$person.time
  TAB<-cbind(time=object$time,hazard=object$hazard,lower=object$low,upper=object$up)
 
#   cat("Smoothing estimates:")
#    cat("\n")        
#          cat(,sep=" ");cat("\n")#Smoothing parameter
#          cat(,sep=" ");cat("\n")#degree of freedom
#          cat(",sep=" ");cat("\n")#overdispersion parameter
  #  print(x$sv2,x$phi,x$df)
  temp<-list(call=object$call,covariate.value=object$cov.value,HazardEstimates=TAB,
             lambda=as.numeric(1/object$sv2),df=as.numeric(object$df),phi=as.numeric(object$phi))
  class(temp) <- "summary.bshazard"
  return(temp)
}


