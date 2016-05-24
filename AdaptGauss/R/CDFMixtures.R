CDFMixtures <- function(Kernels,Means,SDs,Weights = (Means*0)+1,IsLogDistribution = Means*0){
# [CDFGaussMixture,CDFSingleGaussian] = CDFMixtures(Kernels,Means,SDs,Weights,IsLogDistribution);
# gibt die cdf (cumulierte Dichte) einer aus Gauss bzw. LogGauss zusammengesetzten GMM Verteilung zuruck. 
#
# INPUT
# Kernels(1:N)                     at these locations N(Means,SDs)*Weights is used for cdf calcuation
#                                  NOTE:   Kernels are usually (but not necessarily) sorted and unique
# Means(1:L,1)                     Means of Gaussians,  L ==  Number of Gaussians
# SDs(1:L,1)                     estimated Gaussian Kernels = standard deviations
#
# OPTIONAL
# Weights(1:L,1)                   relative number of points in Gaussians (prior probabilities): 
#                                  sum(Weights) ==1, default weight is 1/L
# IsLogDistribution(1:L,1) oder 0  if IsLogDistribution(i)==1, then mixture is lognormal
#                                  default == 0*(1:L)'
#
# OUTPUT
# CDFGaussMixture(1:N,1)           cdf of Sum of SingleGaussians at Kernels
# CDFSingleGaussian(1:N,1:L)       cdf of mixtures at Kernels

# Autor: RG 06/15

#Weights = Weights/sum(Weights)   # Gleichgewichtung
if(length(IsLogDistribution) == 1 && IsLogDistribution == 0) 
  IsLogDistribution = Means*0     # default Gauss


AnzGaussians=length(Means)
CDFSingleGaussian = matrix(0,length(Kernels),AnzGaussians)
CDFGaussMixture =CDFSingleGaussian[,1]

for(g in 1:AnzGaussians){
  if(IsLogDistribution[g] == 1){ # LogNormal 
    sig   = sqrt(log(SDs[g]*SDs[g]/(Means[g]*Means[g])+1))
    mu    = log(abs(Means[g]))-0.5*sig*sig
    if(Means[g] > 0)
    CDFSingleGaussian[,g] = plnorm(Kernels,mu,sig)*Weights[g]
    else 
    CDFSingleGaussian[,g] = (1 -plnorm(-Kernels,mu,sig))*Weights[g]
  }
  else{ # Gaussian
    CDFSingleGaussian[,g] = pnorm(Kernels,Means[g],SDs[g])*Weights[g]    # Gaussian
  }
  CDFGaussMixture =  CDFGaussMixture + CDFSingleGaussian[,g]
}

return(list(CDFGaussMixture = CDFGaussMixture, CDFSingleGaussian = CDFSingleGaussian))
}