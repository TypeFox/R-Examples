bayesPref <-
function(pData=NULL,mcmcL=1000,dirvar=2,calcdic=TRUE,constrain=FALSE,
                    pmpriorLB=1,pmpriorUB=50,ppprior=NULL,dicburn=100,indc=TRUE,
                    pops=TRUE,pminit=NULL,ppinit=NULL,ipinit=NULL,constrainP=NULL,
                    diradd=0.1,univar=2,estip=TRUE,measure="mean"){

  ## Convert pData to matrix
  pData<-as.matrix(pData)
  
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
  
  ## number of choices/approaches for each individual
  nObs<-apply(pData,1,sum)
  
  ## INITIALIZATION
  ## number of individuals
  nInd<-dim(pData)[1]
  ## number of choice catagories
  nCat<-dim(pData)[2]
  ## array for individual prefs
  prefsArray<-array(dim=c(nInd,nCat,mcmcL))
  ## population preference array
  popPrefsL<-vector("list",nPops)
  popAlphasL<-vector("list",nPops)
  popMultL<-vector("list",nPops)
  for (i in 1:nPops){
    popPrefsL[[i]]<-array(dim=c(nCat,mcmcL))
    popAlphasL[[i]]<-array(dim=c(nCat,mcmcL))
    popMultL[[i]]<-numeric(mcmcL)
  }

  ## pr of model given data
  pMD<-numeric(mcmcL)
  ## likelihood, includes conditional prior for dic
  lnlik<-numeric(mcmcL)

  ## set up pop pref prior
  if (is.null(ppprior)==TRUE){
    ppprior<-matrix(rep(1,nCat),nrow=1)
  }
  else if (nCat != length(ppprior)){
    ppprior<-matrix(rep(1,nCat),nrow=1)
    cat("incorrect dimension for ppprior, set to a vector of 1's")
  }
  else{
    ppprior = matrix(ppprior,nrow=1)
  }

  ## initialize values, user supplied or sample from a dirichlet with a similar to observed counts
  if (is.null(ipinit)==FALSE){ # user supplied
    prefsArray[,,1]<-ipinit
  }
  else { # generate default initialization
    for (j in 1:nInd){
      prefsArray[j,,1]<-rdirichlet(1, (pData[j,] + 1))
    }
  }
  ## initialize pop pref, user supplied or sample dirichlet with a similar to mean counts
  if (is.null(ppinit)==FALSE){ # user supplied
    if (is.matrix(ppinit)==FALSE){
      for (i in 1:nPops){	
        popPrefsL[[i]][,1]<-ppinit
      }
    }
    else{
      for (i in 1:nPops){
        popPrefsL[[i]][,1]<-ppinit[i,]
      }
    }	
  }
  else { # generate default initialization
    if (constrain==FALSE){
      for (i in 1:nPops){
        popPrefsL[[i]][,1]<-t(rdirichlet(1,apply(pData,2,mean) + 1))
      }
    }
    else if (constrain==TRUE){
      for (i in 1:nPops){
        popPrefsL[[i]][,1]<-rep((1/nCat),nCat) # all equal
      }
    }
  }
  
  ## uniform for variance
  if (is.null(pminit)==FALSE){ # user supplied
    for (i in 1:nPops){
      popMultL[[i]][1]<-pminit
    }
  }
  else {
    for (i in 1:nPops){
      popMultL[[i]][1]<-runif(1,pmpriorLB,pmpriorUB)
    }
  }
  
  ## population alphas
  for (i in 1:nPops){
    popAlphasL[[i]][,1]<-popPrefsL[[i]][,1] * popMultL[[i]][1]
  }
  
  ## calculate probabilities
  tempprobs<-calcP(step=1,prefsArray=prefsArray,popPrefsL=popPrefsL,
                   popMultL=popMultL,lb=pmpriorLB,ub=pmpriorUB,N=nInd,
                   pprior=ppprior,pData=pData,nPops=nPops,popN=popN,
                   estip=estip)
  lnlik[1]<-tempprobs[1]
  pMD[1]<-tempprobs[2]

  cat("current mcmc step: 1;   p(M|D): ",pMD[1],"\n")
  
  ## BEGIN MCMC LOOP
  for(i in 2:mcmcL){
    ## gibbs sampling update of individual preferences
    if (estip==TRUE){ # estimate individual preference
      for (j in 1:nInd){
        prefsArray[j,,i]<-rdirichlet(1, (pData[j,] + popAlphasL[[popN[j]]][,i-1]))
      }
    }
    ## metropolis-hastings update of population-level preference
    ## proposals, popPrefs by random walk or independence and popMult by uniform
    if (constrain == FALSE){
      if (indc==FALSE){ # random-walk
        for (k in 1:nPops){
          popPrefsL[[k]][,i]<-rdirichlet(1, popPrefsL[[k]][,i-1] * dirvar)
        }
      }
      else { # independence chain
        for (k in 1:nPops){
          popPrefsL[[k]][,i]<-rdirichlet(1, (apply(pData[popN==k,], 2, mean) + diradd) * dirvar)
        }
      }
    }
    else if (constrain == TRUE){ ## pop prefs constant and equal
      for (k in 1:nPops){
        popPrefsL[[k]][,i]<-popPrefsL[[k]][,i-1]
      }
    }

    for (k in 1:nPops){
      popMultL[[k]][i]<-runif(1,popMultL[[k]][i-1] - univar, popMultL[[k]][i-1] + univar)
      if (popMultL[[k]][i] > pmpriorUB | popMultL[[k]][i] < pmpriorLB){
        popMultL[[k]][i]<-popMultL[[k]][i-1]
      }
    }

    for (k in 1:nPops){
      ## popPrefs
      tmpNind<-sum(popN==k)
      tmpprefsArray<-prefsArray[popN==k,,]
      popPrefsL[[k]][,i]<-MHratioPopPref(step=i,nInd=tmpNind,prefsArray=tmpprefsArray,
                                         popPrefs=popPrefsL[[k]],popMult=popMultL[[k]],
                                         popprior=ppprior,dirvar=dirvar,indc=indc,
                                         pData=pData[popN==k,],diradd=diradd,estip=estip)
      ##popMults
      popMultL[[k]][i]<-MratioPopMult(step=i,nInd=tmpNind,prefsArray=tmpprefsArray,
                                      popPrefs=popPrefsL[[k]],popMult=popMultL[[k]],
                                      lb=pmpriorLB,ub=pmpriorUB,estip=estip,pData=pData)

      ## new alphas
      popAlphasL[[k]][,i]<-popPrefsL[[k]][,i] * popMultL[[k]][i]
    }
    ## calc total p. prob
    tempprobs<-calcP(step=i,prefsArray=prefsArray,popPrefsL=popPrefsL,
                     popMultL=popMultL,lb=pmpriorLB,ub=pmpriorUB,N=nInd,
                     pprior=ppprior,pData=pData,nPops=nPops,popN=popN,
                     estip=estip)

    lnlik[i]<-tempprobs[1]
    pMD[i]<-tempprobs[2]
    if (i %% 100 == 0){ ## print every 100 steps
      cat("current mcmc step: ",i,";  p(M|D): ",pMD[i],"\n")
    }
  }

  if (calcdic==TRUE){
    dic<-calcdic(lnL=lnlik,popP=popPrefsL,popM=popMultL,indP=prefsArray,dicburn=dicburn,
                 mcmc=mcmcL,N=nInd,pData=pData,popN=popN,nPops=nPops,estip=estip,measure=measure)
    cat("dic = ",dic,"\n")
  }
  
  ## RETURN RESULTS
  outL<-vector("list",nPops)
  names(outL)<-1:nPops
  for (i in 1:nPops){
    out<-list(prefsArray[popN==i,,],popPrefsL[[i]],popMultL[[i]],lnlik,pMD,dic)
    names(out)<-c("IndPref","PopPref","PopVar","likelihood","pMD","dic")
    if (estip==FALSE){
      out[[1]]<-NULL
    }
    outL[[i]]<-out
  }
  return(outL)
}

