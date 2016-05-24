
.postInit <- function(initialization,Y,Xnorm,Xpois,Xbin,n,k,start.z){
  if (k==1) {
    z <- matrix(rep(1,n),ncol=1)
    return(z)
  }
  if(initialization=="mclust"){
    #z <- gpcm(data=cbind(Y,Xnorm,Xpois,Xbin),G=k,mnames="VVV")$z
    data <- cbind(Y[,1],Xnorm,Xpois,Xbin)
    z <- Mclust(data=data,modelNames=ifelse(ncol(data)==1,"V","VVV"),G=k)$z
  } 
  if(initialization=="kmeans"){
    clusters  <- kmeans(x=cbind(Y,Xnorm,Xpois,Xbin),centers=k)     # clusters on D
    z         <- mclust::unmap(clusters$cluster,1:k)
  } 
  if(initialization=="random.soft"){
    z  <- array(runif(n*k),c(n,k)) # soft posterior probabilities (no-normalized) (n x k) 
    z  <- z/rowSums(z)             # soft posterior probabilities (n x k)
  } 
  if(initialization=="random.hard"){
    z  <- t(rmultinom(n, size = 1, prob=rep(1/k,k)))  # hard posterior probabilities (n x k)
  } 
  if(initialization=="manual"){ # z.start can be both soft and hard initialization
    z  <- start.z      # posterior probabilities (n x k) no-normalized
  }
  z
}