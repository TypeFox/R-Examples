######################################################################
#
# rem.R
#
# Written by Carter T. Butts <buttsc@uci.edu>.
#
# Last Modified 3/08/15
# Licensed under the GNU General Public License version 2 (June, 1991)
#
# Part of the R/relevent package
#
# Contents:
#
#  rem.ord.dev
#  rem.ord.nlp
#  rem.int.dev
#  rem.int.nlp
#  rem
#  print.rem
#  print.summary.rem
#  summary.rem
#
######################################################################

#Calculate deviance (and derivatives) for relational event data, in a manner
#usable by trust.  Arguments are as follows:
#  par - Combined parameter vector
#  pmapl - list of vectors selecting parameters for local stats; ith entry
#    should be the appropriate par element for the ith statistic
#  evl - list of matrices containing observed events; endogenous events are
#    indicated by an integer corresponding to event type in the first column, 
#    while exogenous events are indicated by first-column 0s.  The second
#    column should contain the event timing information, translated so that
#    observation starts at time zero.  For the continuous time case, the last
#    event should be an exogenous event corresponding to the end of the
#    observation period.
#  statsl - list of arrays (event number by event type by statistic) containing
#    statistics for each endogenous event at each event number (i.e., the stats
#    determining the hazard going into the event number in question).  Note that
#    one still has entries here corresponding to exogenous events, as this
#    reflects possible changes in the hazard structure due to the events in
#    question (which will affect the likelihood via the survival functions, at
#    lease in the continuous case).
#  suppl - list of matrices (event number by event type) indicating the support
#    at each event number.  Specifically, suppl[[i]][j,k] should be TRUE iff
#    event type k was possible at event number j in sequence i, otherwise FALSE.
#  ret.deriv - logical; should the derivatives of the deviance (gradient and
#    hessian) also be returned?
#  verbose - logical; should deviance information be printed?
#  Note: this is just the deviance for the data conditional on the local
#    parameters.  For a hierarchical model, this will need to be augmented
#    by information from the higher levels.
rem.ord.dev<-function(par,pmapl,evl,statsl,suppl,ret.deriv,verbose,...){
  np<-length(par)
  val<-0
  if(ret.deriv){
    grad<-rep(0,np)
    hess<-matrix(0,np,np)
  }
  #Calculate contributions from each event set
  for(i in 1:length(statsl)){
    p<-par[pmapl[[i]]]
    np2<-length(p)
    calc<-.C("rem_ord_dev_R",as.double(p),as.integer(np2), as.integer(evl[[i]][,1]),as.integer(NROW(evl[[i]])),as.double(statsl[[i]]),as.integer(dim(statsl[[i]])[2]),as.integer(suppl[[i]]),as.integer(ret.deriv),val=as.double(0),grad=as.double(rep(0.0,np2)),hess=as.double(matrix(0.0,np2,np2)),PACKAGE="relevent",NAOK=TRUE)
    val<-val-2*calc$val
    if(ret.deriv){
      grad[pmapl[[i]]]<-grad[pmapl[[i]]]-2*calc$grad
      hess[pmapl[[i]],pmapl[[i]]]<-hess[pmapl[[i]],pmapl[[i]]]-2*calc$hess
    }
  }
  #Return the result
  if(verbose){
    cat("Printing deviance output\n")
    print(val)
    if(ret.deriv){
      print(grad)
      print(hess)
    }
  }
  if(ret.deriv)
    list(value=val,gradient=grad,hessian=hess)
  else
    val
}


rem.ord.nlp<-function(par,pmapl,evl,statsl,suppl,ppar,ret.deriv,verbose,...){
  #Extract the prior parameters
  mu<-ppar$mu
  sigma<-ppar$sigma
  nu<-ppar$nu
  #Define functions for the derivatives of the log t distribution
  dt<-function(x){           #Why define my own?  No good reason....
    lgamma((nu+1)/2) - (nu+1)/2*log(1+(x-mu)^2/(nu*sigma^2)) - (0.5*(log(nu)+log(pi))+sigma+lgamma(nu/2))
  }
  ddt<-function(x){
    (nu+1)*(mu-x)/(nu*sigma^2+(x-mu)^2)
  }
  dddt<-function(x){
    (nu+1)*((x-mu)^2-nu*sigma^2)/((x-mu)^2+nu*sigma^2)^2
  }
  #First, obtain the deviance
  dev<-rem.ord.dev(par=par,pmapl=pmapl,evl=evl,statsl=statsl,suppl=suppl, ret.deriv=ret.deriv,verbose=verbose,...)
  if(ret.deriv){
    dev$value<-dev$value*0.5               #Convert back to negative likelihood
    dev$gradient<-dev$gradient*0.5
    dev$hessian<-dev$hessian*0.5
    #Now, add the prior
    dev$value<-dev$value-sum(dt(par))
    dev$gradient<-dev$gradient-ddt(par)
    dev$hessian<-dev$hessian-diag(dddt(par))
  }else{
    dev<-dev*0.5                           #Convert to nloglik units
    dev<-dev-sum(dt(par))                  #Add the prior          
  }
  #Return the result
  dev
}


#Interval version of the above -- use timing information, and take exogenous events
#into account directly (via inter-event survival function).
rem.int.dev<-function(par,pmapl,evl,statsl,suppl,ret.deriv,verbose,...){
  np<-length(par)
  val<-0
  if(ret.deriv){
    grad<-rep(0,np)
    hess<-matrix(0,np,np)
  }
  #Calculate contributions from each event set
  for(i in 1:length(statsl)){
    p<-par[pmapl[[i]]]
    np2<-length(p)
    calc<-.C("rem_int_dev_R",as.double(p),as.integer(np2),as.double(evl[[i]]), as.integer(NROW(evl[[i]])),as.double(statsl[[i]]),as.integer(dim(statsl[[i]])[2]),as.integer(suppl[[i]]),as.integer(ret.deriv),val=as.double(0),grad=as.double(rep(0.0,np2)),hess=as.double(matrix(0.0,np2,np2)),PACKAGE="relevent",NAOK=TRUE)
    val<-val-2*calc$val
    if(ret.deriv){
      grad[pmapl[[i]]]<-grad[pmapl[[i]]]-2*calc$grad
      hess[pmapl[[i]],pmapl[[i]]]<-hess[pmapl[[i]],pmapl[[i]]]-2*calc$hess
    }
  }
  #Return the result
  if(verbose){
    cat("Printing deviance output\n")
    print(val)
    if(ret.deriv){
      print(grad)
      print(hess)
    }
  }
  if(ret.deriv)
    list(value=val,gradient=grad,hessian=hess)
  else
    val
}

#Negative log-posterior, currently using a t prior
rem.int.nlp<-function(par,pmapl,evl,statsl,suppl,ppar,ret.deriv,verbose,...){
  #Extract the prior parameters
  mu<-ppar$mu
  sigma<-ppar$sigma
  nu<-ppar$nu
  #Define functions for the derivatives of the log t distribution
  dt<-function(x){           #Why define my own?  No good reason....
    lgamma((nu+1)/2) - (nu+1)/2*log(1+(x-mu)^2/(nu*sigma^2)) - (0.5*(log(nu)+log(pi))+sigma+lgamma(nu/2))
  }
  ddt<-function(x){
    (nu+1)*(mu-x)/(nu*sigma^2+(x-mu)^2)
  }
  dddt<-function(x){
    (nu+1)*((x-mu)^2-nu*sigma^2)/((x-mu)^2+nu*sigma^2)^2
  }
  #First, obtain the deviance
  dev<-rem.int.dev(par=par,pmapl=pmapl,evl=evl,statsl=statsl,suppl=suppl, ret.deriv=ret.deriv,verbose=verbose,...)
  if(ret.deriv){
    dev$value<-dev$value*0.5              #Convert back to negative likelihood
    dev$gradient<-dev$gradient*0.5
    dev$hessian<-dev$hessian*0.5
    #Now, add the prior
    dev$value<-dev$value-sum(dt(par))
    dev$gradient<-dev$gradient-ddt(par)
    dev$hessian<-dev$hessian-diag(dddt(par))
  }else{
    dev<-dev*0.5                          #Convert to negloglik units
    dev<-dev-sum(dt(par))                 #Add the prior
  }
  #Return the result
  dev
}


#Fit a relational event model
rem<-function(eventlist,statslist,supplist=NULL,timing=c("ordinal","interval"),estimator=c("BPM","MLE","BMCMC","BSIR"),prior.param=list(mu=0,sigma=1e3,nu=4),mcmc.draws=1.5e3,mcmc.thin=2.5e1,mcmc.burn=2e3,mcmc.chains=3,mcmc.sd=0.05,mcmc.ind.int=50,mcmc.ind.sd=10,sir.draws=1000,sir.expand=10,sir.nu=4,verbose=FALSE){
  #Multivariate t density (for BSIR)
  dlmvt<-function(x,mu,is,ds,df){
    m<-length(x)
    lgamma((df+m)/2)-m/2*log(pi*df)-lgamma(df/2)-log(abs(ds))/2- (df+m)/2*log(1+((x-mu)%*%is%*%(x-mu))/df)
  }
  #Define an internal convenience function to compute the null deviance in the
  #interval case (constant hazard for all events)
  intNullDev<-function(evl,suppl,nev){
    n<-length(evl)
    stl<-list()
    for(i in 1:n){
      stl[[i]]<-array(suppl[[i]],dim=c(NROW(evl[[i]]),nev,1))
    }
    fitnull<-trust(rem.int.dev,parinit=0,1,100,pmapl= replicate(n,1,simplify=FALSE),evl=evl,statsl=stl,suppl=suppl,ret.deriv=TRUE, verbose=FALSE)
    fitnull$value
  }
  #Make sure that we have lists
  if(!is.list(eventlist))
    eventlist<-list(eventlist)
  if(!is.list(statslist))
    statslist<-list(list(global=statslist))
  #If no timing information given, add it
  for(i in 1:length(eventlist)){
    if((length(eventlist[[i]])>0)&&is.null(dim(eventlist[[i]])))
      eventlist[[i]]<-cbind(eventlist[[i]],1:length(eventlist[[i]]))
  }
  #If no support information given, assume everything is possible all the time
  if(is.null(supplist)){
    supplist<-list()
    for(i in 1:length(eventlist))
      supplist[[i]]<-matrix(TRUE,NROW(eventlist[[i]]), max(dim(statslist[[i]]$global)[2], dim(statslist[[i]]$local)[2]))
  }else if(!is.list(supplist))
    supplist<-list(supplist) 
  #Remove anything with less than two events
  sel<-sapply(eventlist,NROW)>1
  if(any(!sel)){
    warning("Some event lists had less than two events.  Attempting to remove and continue.")
    eventlist<-eventlist[sel]
    statslist<-statslist[sel]
    supplist<-supplist[sel]
  }
  #Initial setup
  listnam<-names(eventlist)                              #Ego names
  if(is.null(listnam))
    listnam<-paste("Ego",1:length(eventlist),sep="")
  pmapl<-list()                                          #Global vars
  statsl<-list()
  if(is.null(statslist[[1]]$global)){
    gp<-0
    np<-0
    parnam<-vector()
  }else{
    gp<-dim(statslist[[1]]$global)[3]
    np<-gp
    parnam<-dimnames(statslist[[1]]$global)[[3]]
    if(is.null(parnam))
      parnam<-paste("t",1:gp,sep="")
  }
  ne<-0
  for(i in 1:length(eventlist)){  #Set up statsl/pmapl and perform sanity checks
    if(is.null(statslist[[i]]$global)){
      gpn<-0
      gpm<-0
      gpe<-0
    }else{
      gpn<-dim(statslist[[i]]$global)[3]
      gpm<-dim(statslist[[i]]$global)[1]
      gpe<-dim(statslist[[i]]$global)[2]
    }
    if(is.null(statslist[[i]]$local)){
      lpn<-0
      lpm<-0
      lpe<-0
    }else{
      lpn<-dim(statslist[[i]]$local)[3]
      lpm<-dim(statslist[[i]]$local)[1]
      lpe<-dim(statslist[[i]]$local)[2]
    }
    epm<-NROW(eventlist[[i]])
    if((epm>1)&&(any(diff(eventlist[[i]][,2])<=0)))
      stop("Simultaneous or improperly ordered events in eventlist ",i,"\n")
    ne<-max(ne,gpe,lpe)
    if(verbose)
      cat("(",i,gp,epm,")",gpn,lpn,gpm,lpm,gpe,lpe,ne,"\n")
    if(gp!=gpn)
      stop("Same number of global statistics must be specified for each event sequence.\n")
    if(((gpn>0)&&(((lpn>0)&&((gpm!=lpm)||(gpe!=lpe))) || (epm!=gpm))) || ((lpn>0)&&(epm!=lpm)))
      stop("Inconsistent event sequence length information for sequence ",i,"\n")
    if(lpn==0){                                #No local stats - use global
      statsl[[i]]<-statslist[[i]]$global
      if(gp>0)
        pmapl[[i]]<-1:gp
      else
        pmapl[[i]]<-NULL
    }else if(gpn==0){                          #No global stats - use local
      statsl[[i]]<-statslist[[i]]$local
      pmapl[[i]]<-np+(1:lpn)
      if(!is.null(dimnames(statslist[[1]]$local)[[3]]))
        parnam<-c(parnam,paste(listnam[i],dimnames(statslist[[1]]$local)[[3]], sep="."))
      else
        parnam<-c(parnam,paste(listnam[i],paste("t", 1:lpn),sep="."))
      np<-np+lpn
    }else{                                      #Both global and local stats
      statsl[[i]]<-array(dim=c(epm,gpe,gpn+lpn))
      statsl[[i]][,,1:gpn]<-statslist[[i]]$global
      statsl[[i]][,,(gpn+1):(gpn+lpn)]<-statslist[[i]]$local
      pmapl[[i]]<-c(1:gp,np+(1:lpn))
      if(!is.null(dimnames(statslist[[1]]$local)[[3]]))
        parnam<-c(parnam,paste(listnam[i],dimnames(statslist[[1]]$local)[[3]], sep="."))
      else
        parnam<-c(parnam,paste(listnam[i],paste("t", 1:lpn),sep="."))
      np<-np+lpn
    }
  }
  #Define the objective function to minimize
  ofun<-switch(paste(match.arg(timing),match.arg(estimator),sep="."),
    "ordinal.MLE"=rem.ord.dev,
    "ordinal.BPM"=rem.ord.nlp,
    "ordinal.BMCMC"=rem.ord.nlp,
    "ordinal.BSIR"=rem.ord.nlp,
    "interval.MLE"=rem.int.dev,
    "interval.BPM"=rem.int.nlp,
    "interval.BMCMC"=rem.int.nlp,
    "interval.BSIR"=rem.int.nlp
  )
  #Fit the model using either posterior simulation or mode finding
  if(verbose)
    cat("Fitting model...\n")
  if(match.arg(estimator)=="BMCMC"){  #Bayesian Posterior Estimation via MCMC
    #Simulate posterior draws using a Metropolis algorithm
    chains<-list()                                    #Chain list
    clen<-ceiling(mcmc.draws/mcmc.chains)             #Chain length
    draws<-matrix(nrow=clen*mcmc.chains,ncol=np)      #Final draws
    geweke<-list()                                    #Gewke diag list
    lpl<-list()                                       #Log posterior list
    lp<-rep(0,clen*mcmc.chains)                       #Final log posteriors
    stime<-proc.time()[3]
    for(i in 1:mcmc.chains){   #Run each chain
      if(verbose)
        cat("\tRunning chain",i,"\n")
      chains[[i]]<-matrix(nrow=clen,ncol=np)
      lpl[[i]]<-rep(0,clen)
      #Initialize with the independence distribution (normal w/fixed sd)
      chains[[i]][1,]<-rnorm(np,sd=mcmc.ind.sd)
      lpl[[i]][1]<--ofun(par=chains[[i]][1,],pmapl=pmapl,evl=eventlist, statsl=statsl,suppl=supplist,ppar=prior.param,ret.deriv=FALSE,verbose=FALSE)
      bc<-1
      dc<-0
      j<-1
      while(j<=clen){
        if(verbose&&((bc+dc)%%1000==0)){
          etime<-proc.time()[3]-stime
          if(bc<mcmc.burn)
            cat("Chain ",i,", Burn-in Iteration ",bc,"/",mcmc.burn,", Elapsed time ",round(etime/60,2)," min (ETA ",round(etime/((i-1)*(mcmc.burn+clen*mcmc.thin)+bc+dc)/60*((mcmc.chains-i)*(mcmc.burn+clen*mcmc.thin)+mcmc.burn-bc+clen*mcmc.thin-dc),2)," min)\n",sep="")
          else
            cat("Chain ",i,", Draw ",j,"/",clen,", Elapsed time ",round(etime/60,2)," min (ETA ",round(etime/((i-1)*(mcmc.burn+clen*mcmc.thin)+bc+dc)/60*((mcmc.chains-i)*(mcmc.burn+clen*mcmc.thin)+mcmc.burn-bc+clen*mcmc.thin-dc),2)," min)\n",sep="")
          cat("\tlp=",round(lpl[[i]][j],4),"; par=(",paste(round(chains[[i]][j,],3),collapse=", "), ")\n",sep="")
        }
        #Create the candidate using a mixture of samplers
        if((bc+dc)%%mcmc.ind.int==0){  #Independence sampler
          cand<-rnorm(np,sd=mcmc.ind.sd)
          ljmprat<-sum(dnorm(chains[[i]][j,],sd=mcmc.ind.sd,log=T)- dnorm(cand,sd=mcmc.ind.sd,log=T))
        }else{                         #Random walk Metropolis sampler
          cand<-chains[[i]][j,]+rnorm(np,sd=mcmc.sd)
          ljmprat<-0 
        }
        cand.lp<--ofun(par=cand,pmapl=pmapl,evl=eventlist,statsl=statsl, suppl=supplist,ppar=prior.param,ret.deriv=FALSE,verbose=FALSE)
        #Decide to keep or reject
        if((lpl[[i]][j]==-Inf)||(cand.lp-lpl[[i]][j]+ljmprat>log(runif(1)))){
          chains[[i]][j,]<-cand                              #Keep
          lpl[[i]][j]<-cand.lp
        }
        #Update the draw counts
        if(bc<=mcmc.burn)
          bc<-bc+1
        else{
          if(dc%%mcmc.thin==0){
            j<-j+1
            if(j<=clen){
              chains[[i]][j,]<-chains[[i]][j-1,]
              lpl[[i]][j]<-lpl[[i]][j-1]
            }
          }
          dc<-dc+1
        }
      }
      draws[((i-1)*clen+1):(i*clen),]<-chains[[i]]
      lp[((i-1)*clen+1):(i*clen)]<-lpl[[i]]
      chains[[i]]<-as.mcmc(chains[[i]])
      geweke[[i]]<-geweke.diag(chains[[i]],frac1=0.15,frac2=0.15)
    }
    if(mcmc.chains>1)
      gelman<-gelman.diag(as.mcmc.list(chains))
    else
      gelman<-NULL
    fit<-list(draws=draws,lp.draws=lp,gelman.rubin=gelman,geweke=geweke, mcmc.burn=mcmc.burn,mcmc.draws=mcmc.chains*clen,mcmc.thin=mcmc.thin,mcmc.chains=mcmc.chains,mcmc.sd=mcmc.sd,mcmc.ind.int=mcmc.ind.int,mcmc.ind.sd=mcmc.ind.sd)
    #Identify the coefficients
    fit$coef<-colMeans(fit$draws)     #Point estimate here is posterior mean
    names(fit$coef)<-parnam
    colnames(fit$draws)<-parnam
  }else{                              #Posterior Mode Estimation/MLE
    #Perform the optimization using trust
    par<-rep(0,np)
    fit<-trust(ofun,parinit=par,1,100,pmapl=pmapl,evl=eventlist, statsl=statsl, suppl=supplist, ppar=prior.param, ret.deriv=TRUE, verbose=verbose)
    #Identify the coefficients
    fit$coef<-fit$argument
    names(fit$coef)<-parnam
  }
  #For BSIR case, need to simulate posterior draws
  if(match.arg(estimator)=="BSIR"){
    if(verbose)
        cat("Taking posterior draws using resampling\n")
    fit$cov.hess<-try(qr.solve(fit$hessian))
    cf<-try(chol(fit$cov.hess))
    if(class(cf)=="try-error")
      stop("Couldn't conduct posterior resampling because estimated covariance matrix was singular; suggest that you examine this model more closely.\n")
    is<-qr.solve(fit$cov.hess)
    ds<-det(fit$cov.hess)
    draws<-matrix(0,sir.draws*sir.expand,np)
    iw<-vector()
    lp<-vector()
    for(i in 1:(sir.draws*sir.expand)){
      draws[i,]<-fit$coef+cf%*%rnorm(np)*sqrt(sir.nu/rchisq(1,sir.nu))
      lp[i]<--ofun(draws[i,],pmapl=pmapl,evl=eventlist,statsl=statsl, suppl=supplist,ppar=prior.param,ret.deriv=FALSE,verbose=FALSE)
      iw[i]<-lp[i]-dlmvt(draws[i,],fit$coef,is,ds,sir.nu)
    }
    iw<-iw-logSum(iw)
    if(any(is.na(iw)|is.nan(iw)|is.na(exp(iw)))||all(exp(iw)==0)){
      warning("Importance weights seem to have some issues.  Carrying on as best we can.\n")
    }
    sel<-sample(1:NROW(draws),sir.draws,replace=TRUE,prob=exp(iw))
    fit$draws<-draws[sel,]
    fit$lp.draws<-lp[sel]
    colnames(fit$draws)<-parnam
    fit$iw<-iw
    fit$coef.mode<-fit$coef
    fit$coef<-colMeans(fit$draws)
    names(fit$coef)<-parnam
  }
  #Final output elaboration and return
  fit$df.model<-np
  fit$df.null<-sum(sapply(eventlist,function(z){sum(z[,1]!=0)}))
  if(match.arg(estimator)=="MLE"){    #If MLE, already have deviance
    fit$residual.deviance<-fit$value
  }else{                              #If Bayesian, have to calculate afresh
    fit$residual.deviance<-switch(match.arg(timing),
      "ordinal"=rem.ord.dev(fit$coef,pmapl=pmapl,evl=eventlist, statsl=statsl, suppl=supplist, ret.deriv=FALSE, verbose=verbose),
      "interval"=rem.int.dev(fit$coef,pmapl=pmapl,evl=eventlist, statsl=statsl, suppl=supplist, ret.deriv=FALSE, verbose=verbose)
    )
  }
  if(match.arg(estimator)%in%c("MLE","BPM")){  #Mode cases
    fit$cov<-try(qr.solve(fit$hessian))
    try(rownames(fit$cov)<-parnam)
    try(colnames(fit$cov)<-parnam)
    if(match.arg(estimator)=="MLE"){
      fit$se<-try(sqrt(diag(fit$cov)))
      try(names(fit$se)<-parnam)
    }else{
      fit$sd<-try(sqrt(diag(fit$cov)))
      try(names(fit$sd)<-parnam)
      fit$logpost<--fit$value
    }
  }else{                                       #MCMC/BSIR cases
    fit$cov<-try(cov(fit$draws))
    try(rownames(fit$cov)<-parnam)
    try(colnames(fit$cov)<-parnam)
    fit$sd<-try(sqrt(diag(fit$cov)))
    try(names(fit$sd)<-parnam)
    fit$logpost<--switch(match.arg(timing),
      "ordinal"=rem.ord.nlp(fit$coef,pmapl=pmapl,evl=eventlist, statsl=statsl, suppl=supplist, ppar=prior.param, ret.deriv=FALSE, verbose=verbose),
      "interval"=rem.int.nlp(fit$coef,pmapl=pmapl,evl=eventlist, statsl=statsl, suppl=supplist, ppar=prior.param, ret.deriv=FALSE, verbose=verbose)
    )
  }
  fit$loglik<--0.5*fit$residual.deviance
  fit$null.deviance<-switch(match.arg(timing),
#    "ordinal"=2*sum((eventlist[,1]>0)*sapply(supplist,function(z){sum(log(rowSums(z)))})),
    "ordinal"=2*sum(sapply(1:length(eventlist),function(z){
      sum((eventlist[[z]][,1]>0)*log(rowSums(supplist[[z]])))
    })),
    "interval"=intNullDev(evl=eventlist,suppl=supplist,nev=ne)
  )
  fit$model.deviance<-fit$null.deviance-fit$residual.deviance
  fit$AIC<-fit$residual.deviance+2*np
  fit$AICC<-fit$AIC+2*np*(np+1)/(fit$df.null-np-1)
  fit$BIC<-fit$residual.deviance+log(fit$df.null)*np
  fit$ordinal<-match.arg(timing)=="ordinal"
  fit$estimator<-match.arg(estimator)
  if(match.arg(estimator)%in%c("BPM","BMCMC","BSIR"))
    fit$prior.param<-prior.param
  class(fit)<-"rem"
  fit
}


#Print method for rem objects
print.rem<-function(x, ...){
    cat("Egocentric Relational Event Model\n")
    cat("\nCoefficients:\n")
    print(x$coef)
    cat("\nNull deviance:",x$null.deviance,"\nResidual deviance:",x$residual.deviance,"\n")
    cat("AIC:",x$AIC,"AICC:",x$AICC,"BIC:",x$BIC,"\n\n")
}


#Print method for summary.rem objects
print.summary.rem<-function (x, ...) 
{
    cat("Egocentric Relational Event Model ")
    if (x$ordinal) 
        cat("(Ordinal Likelihood)\n\n")
    else cat("(Interval Likelihood)\n\n")
    if (is.null(x$cov)) {
        ctab <- matrix(x$coef, ncol = 1)
        rownames(ctab) <- names(x$coef)
        if(x$estimator=="MLE")
          colnames(ctab) <- c("MLE")
        else
          colnames(ctab) <- c("Post.Mode")
        printCoefmat(ctab)
    }
    else {
        ctab <- cbind(x$coef, diag(x$cov)^0.5)
        if(x$estimator%in%c("MLE","BPM")){
          ctab <- cbind(ctab, ctab[, 1]/ctab[, 2])
          ctab <- cbind(ctab, 2 * (1 - pnorm(abs(ctab[, 3]))))
          if(x$estimator=="MLE"){
            colnames(ctab) <- c("MLE", "Std.Err", "Z value", 
              "Pr(>|z|)")
          }else if(x$estimator=="BPM"){
            colnames(ctab) <- c("Post.Mode", "Post.SD", "Z value", 
              "Pr(>|z|)")
          }
        }else{
          ctab <- cbind(ctab, apply(x$draws,2,quantile,0.025), 
            apply(x$draws,2,quantile,0.975),
            apply(x$draws,2,function(z){mean(sign(z)!=sign(mean(z)))}))
          colnames(ctab)<-c("Post.Mean","Post.SD","Q2.5%","Q97.5%","Pr(!sgn(z))")
        }
        rownames(ctab) <- names(x$coef)
        printCoefmat(ctab, P.values = TRUE)
    }
    cat("Null deviance:", x$null.deviance, "on", x$df.null, "degrees of freedom\n")
    cat("Residual deviance:", x$residual.deviance, "on", x$df.null - 
        x$df.model, "degrees of freedom\n")
    cat("\tChi-square:", x$model.deviance, "on", x$df.model, 
        "degrees of freedom, asymptotic p-value", 1 - pchisq(x$model.deviance, 
            x$df.model), "\n")
    cat("AIC:", x$AIC, "AICC:", x$AICC, "BIC:", x$BIC, "\n")
    if(x$estimator%in%c("BPM","BMCMC","BSIR")){
      cat("Log posterior:",x$logpost,"\n")
      cat("Prior parameters:",paste(names(x$prior.param),unlist(x$prior.param),sep="="),"\n")
    }
}


#Summary method for rem objects
summary.rem<-function (object, ...) 
{
    class(object) <- c("summary.rem", class(object))
    object
}

