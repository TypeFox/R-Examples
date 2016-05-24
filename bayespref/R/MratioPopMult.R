MratioPopMult <-
function(step=NULL,nInd=NULL,prefsArray=NULL,popPrefs=NULL,
                        popMult=NULL,lb=NULL,ub=NULL,estip=NULL,pData=NULL){
  ## prior for pop mult
  new_prior<-log(dunif(popMult[step],lb,ub))
  old_prior<-log(dunif(popMult[step-1],lb,ub))

  ## conditional prior prob
  new_conditional<-numeric(nInd)
  old_conditional<-numeric(nInd)
  if (estip==TRUE){
    for (j in 1:nInd){
      new_conditional[j]<-ddirichlet(prefsArray[j,,step], as.vector(popPrefs[,step]) * popMult[step])
      old_conditional[j]<-ddirichlet(prefsArray[j,,step], as.vector(popPrefs[,step]) * popMult[step-1])
    }
  }
  else {
    for (j in 1:nInd){
      new_conditional[j]<-dpolya(pData[j,],as.vector(popPrefs[,step]) * popMult[step])
      old_conditional[j]<-dpolya(pData[j,],as.vector(popPrefs[,step]) * popMult[step-1])
    }
  }
  new_conditional[new_conditional==0]<-10^-323
  old_conditional[old_conditional==0]<-10^-323
  new_conditional[new_conditional==Inf]<-10^308
  old_conditional[old_conditional==Inf]<-10^308
  new_conditional[is.na(new_conditional)==TRUE]<-10^308
  old_conditional[is.na(old_conditional)==TRUE]<-10^308
  newP<-sum(log(new_conditional)) + new_prior
  oldP<-sum(log(old_conditional)) + old_prior

  if (is.na(newP) == TRUE | newP== -Inf){
    newP <- log(10^-323)
    cat("newP:,",newP,"\n",new_conditional,"\n")
  }
  else if (newP== Inf){
    newP<-10^308
  }
  if (is.na(oldP) == TRUE | oldP== -Inf){
    oldP <- log(10^-323)
  }
  
  ## note uniform proposal, so metropolis ratio 
  Mratio<- newP - oldP
  u<-log(runif(1,0,1))
  if (Mratio >= u){
    return(popMult[step])
  }
  else return(popMult[step-1])
}

