single.index.aved<-function(x,y,h=1,kernel="gauss",argd=NULL,take=length(y),seed=1)
{
d<-dim(x)[2]
n<-length(y)

if (take==n) ota<-seq(1,n) else{
   set.seed(seed)
   ota<-ceiling(n*runif(take))
}

if (is.null(argd)){
  grad<-matrix(0,n,d)
  for (ii in 1:d){
     for (jj in ota){
        arg<-x[jj,]
        grad[jj,ii]<-kernesti.der(arg,x,y,h=h,direc=ii,kernel=kernel,vect=FALSE)
     }
  }
  theta<-colMeans(grad)
}
else{
  grad<-matrix(0,d,1)
  for (ii in 1:d){
        grad[ii]<-kernesti.der(argd,x,y,h=h,direc=ii,kernel=kernel,vect=FALSE)
  }
  theta<-grad
}

theta<-theta/sqrt(sum(theta^2))

return(theta)
}

