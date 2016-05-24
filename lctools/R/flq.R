FLQ<-function(Coords,Bandwidth,e,E,Denominator,WType='Bi-square'){
  
  Distances<-dist(Coords)
  Dij <- as.matrix(Distances)
  
  Obs<-length(e)
  
  Wts<-matrix(data=0,nrow=Obs,ncol=Obs)
  
  flq<-matrix(data=0,nrow=Obs,ncol=1)
  lq<-matrix(data=0,nrow=Obs,ncol=1)
  
  for(i in 1:Obs){
    
    w_e<-0.0
    w_E<-0.0
    
    #Get the data and add the distances 
    DataSet<-data.frame(Sector=e,Total=E,DNeighbour=Dij[,i])
    
    #Sort by distance
    DataSetSorted<- DataSet[order(DataSet$DNeighbour),]
    
    #Keep Nearest Neighbours
    SubSet1<-DataSetSorted[1:Bandwidth,]
    
    #Find furthest neighbour
    Kernel_H<-max(SubSet1$DNeighbour)
    
    #Calculate weights
    
    for(j in 1:Bandwidth){
      
      if(WType=='Bi-square'){
        Wts[i,j]<-(1-(SubSet1$DNeighbour[j]/Kernel_H)^2)^2
      }
      else{
        Wts[i,j]<-1
      }
      
      #calculate standardised weights
      Wts[i,]<-Wts[i,]/sum(Wts[i,])
      
      #calculate ratio numerator and denominator separately
      w_e<-w_e+Wts[i,j]*SubSet1$Sector[j]
      
      w_E<-w_E+Wts[i,j]*SubSet1$Total[j]
      
    }
    #calculate FLQ and LQ for each location i
    
    flq[i]<-(w_e/w_E)/Denominator
    lq[i]<-(e[i]/E[i])/Denominator
    
  }
  
  Results<-list(LQ=lq,FLQ=flq)
  
  return(Results)
}