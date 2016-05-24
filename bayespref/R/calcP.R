calcP <-
function(step=NULL,prefsArray=NULL,popPrefsL=NULL,
                popMultL=NULL,lb=NULL,ub=NULL,pprior=NULL,N=NULL,
                pData=NULL,nPops=NULL,popN=NULL,estip=NULL){
  ## hyperprior
  prior<-0
  for (i in 1:nPops){
    prior<- prior + log(dunif(popMultL[[i]][step],lb,ub)) +
      log(ddirichlet(t(popPrefsL[[i]][,step]),pprior))
  }
  ## model with individual preference
  if (estip==TRUE){
    ## conditional prior
    p_cond<-0
    for (i in 1:nPops){
      onepopArray<-prefsArray[popN==i,,]## for one population
      p_cond<-p_cond + sum(log(ddirichlet(onepopArray[,,step], popPrefsL[[i]][,step] *
                                          popMultL[[i]][step])))
    }
    ## L(D|conditional prior)
    L<-0
    for(i in 1:N){ ## returns log, hence addition
      L<- L + dmultinom(pData[i,],size=sum(pData[i,]),prob=prefsArray[i,,step],
                        log=TRUE)
    }

    if (is.na(p_cond) == TRUE | p_cond == -Inf){
      p_cond <-  log(10^-323)
    }
    else if (p_cond == Inf){
      p_cond <- log(10^308)
    }
    if (is.na(L) == TRUE | L == -Inf){
      L <-  log(10^-323)
    }
    else if (L == Inf){
      L<- log(10^308)
    }

    
    lnL<- p_cond + L
    pMD<- p_cond + L + prior
  }

  ## model without individual preference
  else{
    L<-0
    for(i in 1:N){
      curpop<-popN[i]
      L<- L + log(dpolya(pData[i,],alpha=popPrefsL[[curpop]][,step] * popMultL[[curpop]][step]))
    }
    if (is.na(L) == TRUE | L == -Inf){
      L <-  log(10^-323)
    }
    lnL<- L
    pMD<- L + prior
  }
  res<-c(lnL,pMD)
  return(res)
}

