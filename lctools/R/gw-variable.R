gw_variable<-function(Coords,InputVariable){
  
  Distances<-dist(Coords)
  m <- as.matrix(Distances)
  
  Obs<-length(InputVariable)
  
  Regional<-matrix(data=NA,nrow=Obs,ncol=1)
  a<-matrix(data=NA,nrow=Obs,ncol=1)
  p<-matrix(data=NA,nrow=Obs,ncol=1)
  
  for(i in 1:Obs){
    for(j in 1:Obs){
      if(j!=i){
        a[i]<-(InputVariable[j]/InputVariable[i])*(1/(m[i,j])^2)
        p[i]<-(1/(m[i,j])^2)
      }
    }
    Regional[i]<-sum(a[i])/sum(p[i])
  }
  return(Regional)
}