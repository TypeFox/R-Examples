#----------------------------------------------------------
#
# Calculate CRPS decomposition.
# from alphas, betas and heavisides
# This function is called by crpsDecomposition
# Returns:   
#         CRPS: mean CRPS
#         Reli: The reliability term of the CRPS
#         CRPSpot: The potential CRPS    
# Author: Ronald Frenette, Severe Weather Lab, Quebec region
#         Jun 2009
# 
#-----------------------------------------------------------
crpsFromAlphaBeta<-function(alpha,beta,heaviside0,heavisideN)
{
  nMember=dim(alpha)[2] -1  
  Reli<-0
  CRPSpot<-0
  for (i in 0:nMember)
  {
    index<-i+1

    meanoi<-0
    meangi<-0

    #Outlier
    if (i==0)
    {
      meanbeta<-mean(beta[,index])
      meanoi<-mean(heaviside0)
      if (meanoi != 0)
        meangi<-meanbeta/meanoi
    }
    if (i==nMember)
    {
      meanoi<-mean(heavisideN)
      meanalpha<-mean(alpha[,index])
      if (meanoi != 1)
        meangi<-meanalpha/(1-meanoi)
    }

    #Non outliers
    if (i>0 & i<nMember)
    {
      meanbeta<-mean(beta[,index])
      meanalpha<-mean(alpha[,index])
      meanoi<-meanbeta/(meanalpha+meanbeta)
      meangi<-meanalpha+meanbeta
    }


    pi<-i/nMember

    Reli<- Reli + meangi * (meanoi-pi) * (meanoi-pi)
    CRPSpot<- CRPSpot + meangi * meanoi * (1.0 - meanoi)
  }

  CRPS<-Reli+CRPSpot
  return(list(CRPS=CRPS,CRPSpot=CRPSpot,Reli=Reli))
}



#----------------------------------------------------------
#
# Calculate CRPS decomposition.
# from observations and ensemble forecast
# Returns:   
#         CRPS: mean CRPS
#         Reli: The reliability term of the CRPS
#         CRPSpot: The potential CRPS   
#         alpha: vector of alpha 
#         beta: vector of beta 
#         heaviside0 vector of heaviside of first outliers
#         heaviside0 vector of heaviside of last outliers
# Author: Ronald Frenette, Severe Weather Lab, Quebec region
#         Jun 2009
# 
#-----------------------------------------------------------
crpsDecomposition<-function(obs,eps)
{
  nMember=dim(eps)[2]
  nObs<-length(obs)
  alpha<-rep(0,nObs*(nMember+1))
  beta<-rep(0,nObs*(nMember+1))
  heaviside0<-rep(0,nObs)
  heavisideN<-rep(0,nObs)
  dim(alpha)<-c(nObs,nMember+1)
  dim(beta)<-c(nObs,nMember+1)
  prev<-t(apply(eps,1,sort))

  # Calculate alpha and beta of observation
  # heaviside for the two outliers

  #1) Beta and alpha for Outliers
  index<-which(obs < prev[,1])
  beta[index,1]<-prev[index,1]-obs[index]
  index<-which(obs > prev[,nMember])
  alpha[index,nMember+1]<-obs[index]-prev[index,nMember]

  #2) Heavisides for Outliers
  index<-which(obs <= prev[,1])
  heaviside0[index]<-1
  index<-which(obs <= prev[,nMember])
  heavisideN[index]<-1

  #3) Non outlier

  for (i in 1:(nMember-1))
  {
    index<-which(obs > prev[,i+1])
    alpha[index,i+1]<-prev[index,i+1]-prev[index,i]
    index<-which(obs < prev[,i])
    beta[index,i+1]<-prev[index,i+1]-prev[index,i]
    index<-which((prev[,i+1] > obs) & (obs > prev[,i]))
    alpha[index,i+1]<-obs[index]-prev[index,i]
    beta[index,i+1]<-prev[index,i+1]-obs[index]
  }

  crps<-crpsFromAlphaBeta(alpha,beta,heaviside0,heavisideN)
  return(list(CRPS=crps$CRPS,CRPSpot=crps$CRPSpot,Reli=crps$Reli,alpha=alpha,beta=beta,heaviside0=heaviside0,heavisideN=heavisideN))
}
