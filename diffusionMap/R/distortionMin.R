distortionMin <- function(X,phi0,K,c0,epsilon=0.001){

  n=dim(X)[1]
  c=c0
  oldD=Inf

  MaxIter=1000

  for(ii in 1:MaxIter){ #K-means loop
    DX=c()
    for(jj in 1:K){
      dX = X-matrix(1,n,1)%*%c[jj,] # n by p
      DX=cbind(DX,apply(dX*dX,1,sum))# n by K
    }
    S = apply(DX,1,which.min) # new labels
    Dtmp = apply(DX,1,min) # new dists. to centroids

    for(jj in 1:K){ # check for empty clusters
      ind=which(S==jj)
      if(length(ind)==0){# if cluster jj empty
        S[which.max(Dtmp)]=jj # give it the pt. furthest from its centroid
        Dtmp[which.max(Dtmp)]=0 # and set its distance to 0
      }
    }

    for(jj in 1:K){# update diffusion centroids
      ind=which(S==jj)
      c[jj,]=t(phi0[ind])%*%X[ind,]/sum(phi0[ind]) #centroid of cluster jj
    }
    D = Dtmp%*%phi0 # distortion

    if((oldD-D)/D < epsilon){ #stopping criterion
      break
    }
    oldD=D
  }
  if(ii==MaxIter) print('Maximum # of iterations reached')
  
  return(list(S=S,c=c,D=D,DX=DX))
}
