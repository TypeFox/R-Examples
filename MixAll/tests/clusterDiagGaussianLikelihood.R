# TODO: Add comment
#
# Author: iovleff
#-----------------------------------------------------------------------
library(MixAll)
data(iris)

gauss_model <- clusterDiagGaussian( iris[1:4], nbCluster = 3, modelNames= c("gaussian_pk_sjk")
                                  , strategy = clusterFastStrategy())

data<-gauss_model@component@data
nbCluster <- gauss_model@nbCluster
prop <- gauss_model@pk
mean <- gauss_model@component@mean
sigma <- gauss_model@component@sigma

nbSample <- nrow(data)
nbVariable <- ncol(data)

f <-vector(length=nbSample)
lnComp <- vector(length=nbCluster)

for (i in 1:nbSample)
{
  for (k in 1:nbCluster)
  { lnComp[k] = log(prop[k]) + sum(dnorm(data[i,], mean[k,], sigma[k,],log=TRUE)); }
  lmax <- max(lnComp)

  lnComp =  lnComp -lmax;

  f[i] = log(sum(exp(lnComp))) + lmax;
}

if( abs(sum(f) - gauss_model@lnLikelihood) < 1.e-15 )
{ print ("clusterDiagGaussianLikelihood failed")
} else
{ print ("clusterDiagGaussianLikelihood successful") }
