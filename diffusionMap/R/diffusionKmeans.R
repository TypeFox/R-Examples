diffusionKmeans <- function(dmap, K, params=c(), Niter=10, epsilon=0.001){

  n=dim(dmap$X)[1]
  D=Inf # max. distortion

  for(ii in 1:Niter){# run Niter K-means loops
    print(paste('Iteration',ii,'of',Niter))
    c0 = dmap$X[sample(1:n,K),] # choose random initial centroids
    kmeans = distortionMin(dmap$X,dmap$phi0,K,c0,epsilon) # K-means loop
    if(kmeans$D<D){ # keep best result
      D = kmeans$D
      DX = kmeans$DX
      part = kmeans$S
      cent = kmeans$c
    }
  }

  if(!is.null(params)){# if parameters are given for each data point, compute the centroid parameters
    npar = dim(params)[2]
    npar = ifelse(is.null(npar),1,npar)
    params = matrix(params,n,npar)
    centparams=matrix(0,K,npar)
    for(jj in 1:K){
      ind=which(part==jj)
      centparams[jj,] = t(dmap$phi0[ind])%*%params[ind,]/sum(dmap$phi0[ind])
    }
    return(list(part=part,cent=cent,D=D,DX=DX,centparams=centparams))
  }
  else{
    return(list(part=part,cent=cent,D=D,DX=DX))
  }
}
