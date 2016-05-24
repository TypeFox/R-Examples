l.moransI<-function(Coords, Bandwidth, x, WType='Binary'){

  Distances<-dist(Coords)
  Dij <- as.matrix(Distances)
  
  #Observations
  Obs<-length(x)
  
  if(Bandwidth>=Obs){
    Bandwidth<-Obs-1
    msg<-cat("Bandwidth set to:",Bandwidth)}
    
  #Mean(x)
  mean.x=mean(x)
  sd.x=sd(x)
   
  #Inference
  n<-Obs
  m4<-sum((x - mean.x)^4)/n
  m2<-sum((x - mean.x)^2)/n
  b2<-m4/(m2^2)
 
  #Local moran
  l.moran<-matrix(data=NA,nrow=Obs,ncol=8)
  l.moran_nom<-matrix(data=0,nrow=Obs,ncol=1)
  moran_denom<-0.0
  
  Wts<-matrix(data=0,nrow=Obs,ncol=Obs)
  
  for(i in 1:Obs){    
    
    l.moran[i,6]<-(x[i] - mean.x)/sd.x
    wXj<-0.0
    
    #Get the data and add the distances 
    DataSet<-data.frame(x,DNeighbour=Dij[,i])
    
    #Sort by distance
    DataSetSorted<- DataSet[order(DataSet$DNeighbour),]
    
    #Keep Nearest Neighvbours
    SubSet1<-DataSetSorted[2:(Bandwidth+1),]
    
    #Find furthest neighbour
    Kernel_H<-max(SubSet1$DNeighbour)
    
    #Calculate weights
    for(j in 1:Obs){
     
      if (DataSet$DNeighbour[j] > Kernel_H){
        Wts[i,j]<-0 
        }
      else{
        if(WType=='Bi-square'){
          Wts[i,j]<-(1-(DataSet$DNeighbour[j]/Kernel_H)^2)^2}
        else{
          Wts[i,j]<-1/Bandwidth} #/Bandwidth
        }
      
      if (j!=i){
        l.moran_nom[i,1]<-l.moran_nom[i,1] + (Wts[i,j]*(x[j]-mean.x))
        wXj<-wXj + (Wts[i,j]*x[j])
      }
       else{
        Wts[i,j]<-0}
       }
    l.moran[i,1]<-(x[i]-mean.x)*l.moran_nom[i,1]
    l.moran[i,7]<-wXj
    
    moran_denom<-moran_denom + ((x[i]-mean.x)*(x[i]-mean.x))
    
    #Inference
   
    #local expected I E(I)
    l.moran[i,2]<--sum(Wts[i,])/(n-1)
    
    #local Var(I)
    a1<-sum(Wts[i,]*Wts[i,])
    A1<-(n-b2)/(n-1)
    a2<-sum(sum(crossprod(Wts[i,])))
    A2<-(2*b2-n)/((n-1)*(n-2))
    
    l.moran[i,3]<-a1*A1 + a2*A2 - ((sum(Wts[i,])^2)/((n-1)^2))
    
   }
  
  Denom <-moran_denom / Obs
  
  #local I
  l.moran[,1]<-l.moran[,1]/Denom
  
  #local Z
  l.moran[,4]<-(l.moran[,1]-l.moran[,2])/sqrt(l.moran[,3])
  
  #local p
  l.moran[,5]<-2*pnorm(-abs(l.moran[,4]))
  
  mean.wXj<-mean(l.moran[,7])
  sd.wXj<-sd(l.moran[,7])
  
  l.moran[,7]<-(l.moran[,7]-mean.wXj)/ sd.wXj
  
  for(i in 1:Obs){
    if (l.moran[i,5]>0.05) {l.moran[i,8]<-0}
    else {
      if(l.moran[i,6]>0 && l.moran[i,7]>0) {l.moran[i,8]<-1}
      else{
        if(l.moran[i,6]<0 && l.moran[i,7]<0) {l.moran[i,8]<-2}
        else{
          if(l.moran[i,6]<0 && l.moran[i,7]>0) {l.moran[i,8]<-3}
          else{
            if(l.moran[i,6]>0 && l.moran[i,7]<0) {l.moran[i,8]<-4}
          }}}
    }
  }
 
  Results<-data.frame(ID=c(1:Obs),Ii=l.moran[,1], Ei=l.moran[,2],Vi=l.moran[,3], Zi=l.moran[,4],p.value=l.moran[,5], Xi=l.moran[,6],wXj=l.moran[,7], Cluster=l.moran[,8])
  
  return(Results)
}