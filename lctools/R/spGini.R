spGini<-function(Coords,Bandwidth,x,WType='Binary'){

  Dij <- as.matrix(dist(Coords))
  
  Obs<-length(x)
  
  gGini_nom<-0.0
  gwGini_nom<-0.0
  nsGini_nom<-0.0
   
  Wts<-matrix(data=0,nrow=Obs,ncol=Obs)
  
  for(m in 1:Obs){
    #Get the data and add the distances 
    DataSet<-data.frame(x,DNeighbour=Dij[,m])
    
    #Sort by distance
    DataSetSorted<- DataSet[order(DataSet$DNeighbour),]
    
    #Keep Nearest Neighvbours
    SubSet1<-DataSetSorted[1:Bandwidth,]
    
    #Find furthest neighbour
    Kernel_H<-max(SubSet1$DNeighbour)
    
    #Calculate weights
    
    for(j in 1:Obs){
     
      if (DataSet$DNeighbour[j] > Kernel_H){
        Wts[m,j]<-0 
        }
      else{
        if(WType=='Bi-square'){
          Wts[m,j]<-(1-(DataSet$DNeighbour[j]/Kernel_H)^2)^2}
        else{
          Wts[m,j]<-1}
      }
      
      if (j==m){
        Wts[m,j]<-0
      }
    }
    if(WType=='RSBi-square'){Wts[m,]<-Wts[m,]/sum(Wts[m,])}
    
    for(j in 1:Obs){
      
    gGini_nom<-gGini_nom+abs(x[m]-x[j])
    
    gwGini_nom<-gwGini_nom+Wts[m,j]*abs(x[m]-x[j])
    
    nsGini_nom<-nsGini_nom+(1-Wts[m,j])*abs(x[m]-x[j])
    }
    
  }
  
  Denom<-2*(Obs^2)*mean(x)
  
  Gini<-gGini_nom/Denom
  gwGini=gwGini_nom/Denom
  nsGini=(gGini_nom - gwGini_nom)/Denom
  
  return(c(Gini=Gini,gwGini=gwGini,nsGini=nsGini,gwGini.frac=gwGini/Gini, nsGini.frac=nsGini/Gini))
  
}