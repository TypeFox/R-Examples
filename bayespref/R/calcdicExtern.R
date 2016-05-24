calcdicExtern <-
function(pData=NULL,prefres=NULL,pops=TRUE,dicburn=0,estip=TRUE,
                        measure="median",constrainP=NULL){
  lnL<-prefres[[1]]$likelihood
  dicvec<--2 * lnL
  if (measure=="mean"){
    dbar<-mean(dicvec)
  }
  else if (measure=="median"){
    dbar<-median(dicvec)
  }

  mcmc<-length(prefres[[1]]$pMD)

  ## DETERMINE WHETHER THERE ARE MULTIPLE POPULATIONS
  ## grab population information if it was supplied, otherwise assume one population
  if (pops==TRUE){
    popN<-pData[,1] # vector of population numbers for each individual
    pData<-pData[,-1] # pData is now missing population information
    if (is.null(constrainP)==FALSE){ # impose constraints for shared parameters
      for (i in 1:max(popN)){
        popN[popN==i]<-constrainP[i]
      }
    }
    nPops<-max(popN) # number of populations being modeled
  }
  else { # if pops is false assume one population
    nPops<-1
    popN<-rep(1,dim(pData)[1])
  }
  
  if(estip==TRUE){ ## model including individual preferences
    p_cond<-0
    L<-0
    if (measure=="mean"){
      for (k in 1:nPops){ ## calculate conditional prob of ind preferences for each
        ## population based on mean or median from posterior distribution for parameters
        pp<-apply(prefres[[k]]$PopPref[,dicburn:mcmc],1,mean)
        pm<-mean(prefres[[k]]$PopVar[dicburn:mcmc])
        alpha<-pp * pm
        ip<-apply(prefres[[k]]$IndPref[,,dicburn:mcmc],c(1,2),mean)
        ip<-ip/apply(ip,1,sum)
        p_cond<-p_cond + sum(log(ddirichlet(ip,alpha)))
        onepData<-pData[popN==k,]
        N<-dim(onepData)[1]
        for(i in 1:N){ ## returns log, hence addition
          L<- L + dmultinom(onepData[i,],size=sum(onepData[i,]),prob=ip[i,],
                            log=TRUE)
        }
      }
    }
    else if (measure=="median"){
      for (k in 1:nPops){ ## calculate conditional prob of ind preferences for each
        ## population based on mean or median from posterior distribution for parameters
        pp<-apply(prefres[[k]]$PopPref[,dicburn:mcmc],1,median)
        pm<-median(prefres[[k]]$PopVar[dicburn:mcmc])
        alpha<-pp * pm
        ip<-apply(prefres[[k]]$IndPref[,,dicburn:mcmc],c(1,2),median)
        ip<-ip/apply(ip,1,sum)
        p_cond<-p_cond + sum(log(ddirichlet(ip,alpha)))
        onepData<-pData[popN==k,]
        N<-dim(onepData)[1]
        for(i in 1:N){ ## returns log, hence addition
          L<- L + dmultinom(onepData[i,],size=sum(onepData[i,]),prob=ip[i,],
                            log=TRUE)
        }
      }
    }
    dhat<- -2 * (p_cond + L)
  }

  else { ## model excluding individual preferences
    L<-0
    N<-dim(pData)[1]
    if (measure=="mean"){
      for(i in 1:N){
        curpop<-popN[i]
        pp<-apply(prefres[[curpop]]$PopPref[,dicburn:mcmc],1,mean)
        pm<-mean(prefres[[curpop]]$PopVar[dicburn:mcmc])
        alpha<-pp * pm
        L<- L + log(dpolya(pData[i,],alpha=alpha))
      }
    }
    else if (measure=="median"){
      for(i in 1:N){
        curpop<-popN[i]
        pp<-apply(prefres[[curpop]]$PopPref[,dicburn:mcmc],1,median)
        pm<-median(prefres[[curpop]]$PopVar[dicburn:mcmc])
        alpha<-pp * pm
        L<- L + log(dpolya(pData[i,],alpha=alpha))
      }
    }
    dhat<- -2 * L
  }

  
  pD<-dbar - dhat
  dic<-pD + dbar
  return(dic)
}

