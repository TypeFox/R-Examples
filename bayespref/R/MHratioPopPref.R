MHratioPopPref <-
function(step=NULL,nInd=NULL,prefsArray=NULL,popPrefs=NULL,
                         popMult=NULL,popprior=NULL,dirvar=NULL,indc=NULL,
                         pData=NULL,diradd=NULL,estip=NULL){
  
  ## prior for pop pref
  new_prior<-log(ddirichlet(popPrefs[,step],popprior))
  old_prior<-log(ddirichlet(popPrefs[,step-1],popprior))
  
  ## conditional prior prob
  new_conditional<-numeric(nInd)
  old_conditional<-numeric(nInd)
  if (estip==TRUE){
    for (j in 1:nInd){
      new_conditional[j]<-ddirichlet(prefsArray[j,,step], as.vector(popPrefs[,step]) * popMult[step-1])
      old_conditional[j]<-ddirichlet(prefsArray[j,,step], as.vector(popPrefs[,step-1]) * popMult[step-1])
    }
  }
  else {
    for (j in 1:nInd){
      new_conditional[j]<-dpolya(pData[j,],as.vector(popPrefs[,step]) * popMult[step-1])
      old_conditional[j]<-dpolya(pData[j,],as.vector(popPrefs[,step - 1]) * popMult[step-1])
    }
  }
  ## clean up probilities replacing
  new_conditional[new_conditional==0]<-10^-323
  old_conditional[old_conditional==0]<-10^-323
  new_conditional[new_conditional==Inf]<-10^308
  old_conditional[old_conditional==Inf]<-10^308
  new_conditional[is.na(new_conditional)==TRUE]<-10^308
  old_conditional[is.na(old_conditional)==TRUE]<-10^308

  ## conditional prior and hyperprior
  newP<-sum(log(new_conditional)) + new_prior
  oldP<-sum(log(old_conditional)) + old_prior
  if (is.na(newP) == TRUE | newP== -Inf){
    cat("newP:,",newP,"\n",new_conditional,"\n",prefsArray[,1,step],"\n",prefsArray[,2,step],"\n\n")
    newP <- log(10^-323)
  }
  else if (newP== Inf){
    newP<-10^308
  }
  if (is.na(oldP)==TRUE | oldP== -Inf){
    cat("oldP:,",oldP,"\n")
    oldP <- log(10^-323)
  }
  
  ## proposal probs.
  if (indc==FALSE){ # random-walk
    old_new<-log(ddirichlet(popPrefs[,step-1],popPrefs[,step] * dirvar))
    new_old<-log(ddirichlet(popPrefs[,step],popPrefs[,step - 1] * dirvar))
  }
  else { # independence chain
    old_new<-log(ddirichlet(popPrefs[,step-1],(apply(pData, 2, mean) + diradd) * dirvar))
    new_old<-log(ddirichlet(popPrefs[,step], (apply(pData, 2, mean) + diradd) * dirvar))
  }
  MHratio<- (newP - new_old) - (oldP - old_new)
  u<-log(runif(1,0,1))
  if (MHratio >= u){
    return(popPrefs[,step])
  }
  else return(popPrefs[,step-1])
}

