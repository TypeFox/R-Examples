# TODO: Add comment
#
# Author: iovleff
#-----------------------------------------------------------------------
library(MixAll)
data(iris)

gamma_model <- clusterGamma( iris[1:4], nbCluster = 3, modelNames = c("gamma_pk_ajk_bjk")
                           , strategy = clusterFastStrategy())

data<-gamma_model@component@data
nbCluster <- gamma_model@nbCluster
prop <- gamma_model@pk
shape <- gamma_model@component@shape
scale <- gamma_model@component@scale

nbSample <- nrow(data)
nbVariable <- ncol(data)

f <-vector(length=nbSample)
lnComp <- vector(length=nbCluster)

for (i in 1:nbSample)
{
  for (k in 1:nbCluster)
  { lnComp[k] = log(prop[k]) + sum(dgamma(data[i,], shape=shape[k,], scale=scale[k,],log=TRUE)); }
  lmax <- max(lnComp)

  for (k in 1:nbCluster)
  { lnComp[k] =  lnComp[k] -lmax;}

  f[i] = log(sum(exp(lnComp))) + lmax;
}
if (abs(sum(f) - gamma_model@lnLikelihood) < 1.e-15)
{ print ("clusterGammaLikelihood failed")
} else{ print ("clusterGammaLikelihood successful")  }
