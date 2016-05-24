
dhist<-function(x, z, pz){
  fx<-rep(0,length(x))
  
  h <- rep(0, length(z) -1)
  for(i in 1:length(h)){
    h[i]<-(pz[i+1] - pz[i]) / (z[i+1]-z[i])
  }
  
  nz<-length(z)
  
  for(i in 1:length(x)){
    index<- (x[i]<=z[2:nz]) & (x[i]>z[1:(nz-1)])
    if(sum(index)>0){
      fx[i] <- h[index]
    }
  }
  fx
}



phist<-function(x,z,pz){
  z<-c(-exp(100), z, exp(100))
  pz<-c(0,pz,1)
  approx(z,pz,x)$y
}

qhist<-function(q,z,pz){
  approx(pz,z,q)$y
}