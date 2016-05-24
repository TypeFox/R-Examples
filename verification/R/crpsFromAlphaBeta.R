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
      if ((meanalpha + meanbeta) != 0) {

         meanoi <- meanbeta/(meanalpha + meanbeta)

       }
      # meanoi<-meanbeta/(meanalpha+meanbeta)
      meangi<-meanalpha+meanbeta
    }


    pi<-i/nMember

    Reli<- Reli + meangi * (meanoi-pi) * (meanoi-pi)
    CRPSpot<- CRPSpot + meangi * meanoi * (1.0 - meanoi)
  }

  CRPS<-Reli+CRPSpot
  return(list(CRPS=CRPS,CRPSpot=CRPSpot,Reli=Reli))
}

