#this function makes simulation data set for single stage
make_classification<-function(n_cluster,pinfo,pnoise,n_sample,centroids=0){
  mu=numeric(pinfo)
  Sigma=diag(pinfo)
  if(is.vector(centroids)){  
    centroids=5*mvrnorm(n_cluster,mu,Sigma)}

X=matrix(0,n_sample,pinfo+pnoise)
y=numeric(n_sample)
z=numeric(n_sample)

stopp=1
A=2*rbinom(n_sample,1,0.5)-1


for (k in 1:n_cluster){
start=stopp
stopp=stopp+n_sample/n_cluster
y[start:(stopp-1)]=2*(k%%2)-1
X[start:(stopp-1),1:pinfo]=matrix(1,n_sample/n_cluster,1)%*%centroids[k,]+mvrnorm(n_sample/n_cluster,mu,Sigma)
z[start:(stopp-1)]=1.5*(2*(k%%2)-1)*A[start:(stopp-1)]+rnorm(n_sample/n_cluster)
}

mun=numeric(pnoise)
Sigman=diag(pnoise)
X[,pinfo+(1:pnoise)]=mvrnorm(n_sample,mun,Sigman)
example=list(X=X,A=A,y=y,R=z,centroids=centroids)
class(example)<-'examplesingle'
return(example)
}



make_2classification<-function(n_cluster,pinfo,pnoise,n_sample,centroids=0){
mu=numeric(pinfo)
Sigma=diag(pinfo)
if(is.vector(centroids)){  
  centroids=5*mvrnorm(n_cluster,mu,Sigma)}

X=matrix(0,n_sample,pinfo+pnoise)
y=list()
z=list()
A=list()
for (i in 1:2){
y[[i]]=rep(1,n_sample)
z[[i]]=rep(1,n_sample)
}

stopp=1
A[[1]] = 2*rbinom(n_sample,1,0.5)-1
A[[2]] = 2*rbinom(n_sample,1,0.5)-1



for (k in 1:n_cluster){
  start=stopp
  stopp=stopp+n_sample/n_cluster
  y[[1]][start:(stopp-1)]=2*(k%%2)-1
  y[[2]][start:(stopp-1)]=2*(floor(k/2)%%2)-1
  X[start:(stopp-1),1:pinfo]=matrix(1,n_sample/n_cluster,1)%*%centroids[k,]+mvrnorm(n_sample/n_cluster,mu,Sigma)
  z[[2]][start:(stopp-1)]=y[[1]][start:(stopp-1)]*A[[1]][start:(stopp-1)]+y[[2]][start:(stopp-1)]*A[[2]][start:(stopp-1)]+rnorm(n_sample/n_cluster)
}
z[[1]]=rep(0,n_sample)
mun=numeric(pnoise)
Sigman=diag(pnoise)
X[,pinfo+(1:pnoise)]=mvrnorm(n_sample,mun,Sigman)
example=list(X=X,A=A,y=y,R=z,centroids=centroids)
class(example)<-'example2'
return(example)
}


