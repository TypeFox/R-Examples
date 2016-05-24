moransI<-function(Coords,Bandwidth,x,WType='Binary'){

  Distances<-dist(Coords)
  Dij <- as.matrix(Distances)
  
  Obs<-length(x)
  
  if(Bandwidth>=Obs){
    Bandwidth<-Obs-1
    msg<-cat("Bandwidth set to:",Bandwidth)}
  
  moran_nom<-0.0
  moran_denom<-0.0
  mean.x=mean(x)
   
  Wts<-matrix(data=0,nrow=Obs,ncol=Obs)
  
  for(i in 1:Obs){
    #Get the data and add the distances 
    DataSet<-data.frame(x,DNeighbour=Dij[,i])
    
    #Sort by distance
    DataSetSorted<- DataSet[order(DataSet$DNeighbour),]
    
    #Keep Nearest Neighbours
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
          Wts[i,j]<-1}
        }
      
      if (j!=i){
         moran_nom<-moran_nom + Wts[i,j]*(x[i]-mean.x)*(x[j]-mean.x)
         }
      }
    moran_denom<-moran_denom + ((x[i]-mean.x)*(x[i]-mean.x))
  }
  
  diag(Wts)<-0
  sum.w=sum(Wts)
  
  Nom<-Obs*moran_nom
  
  Denom<-sum.w*moran_denom
  
  morans.I<-Nom/Denom
  
  #Inference
  #Observations
  n<-Obs
  
  #Expected I E(I)
  E.I<-(-1)/(Obs-1)
  
    
  S0<-sum.w
  
  S1<-(1/2.0) * sum((Wts + t(Wts))^2)
  
  S2<-sum((apply(Wts, 1, sum) + apply(Wts, 2, sum))^2)
  
  b2<-(sum((x - mean.x)^4)/n)/((sum((x - mean.x)^2)/n)^2)
    
  #Var(I)
  Var.I.resampling<-(((n^2) * S1 - n*S2 + 3 * (S0^2))/(((n^2)-1)*(S0^2)))-(E.I^2)
  Var.I.randomization<-(n*((n^2-3*n+3)*S1-n*S2+3*S0^2))/((n-1)*(n-2)*(n-3)*S0^2)-(b2*((n^2-n)*S1-2*n*S2+6*S0^2))/((n-1)*(n-2)*(n-3)*S0^2)-(E.I^2)
  
  Z.I.resampling<-(morans.I-E.I)/sqrt(Var.I.resampling)
  Z.I.randomization<-(morans.I-E.I)/sqrt(Var.I.randomization)

  pv.resampling <- 2*pnorm(-abs(Z.I.resampling))
  pv.randomization <- 2*pnorm(-abs(Z.I.randomization))
  
  Results<-list(W=Wts,Morans.I=morans.I, Expected.I=E.I, z.resampling=Z.I.resampling,z.randomization=Z.I.randomization, p.value.resampling=pv.resampling,p.value.randomization=pv.randomization)
  
  return(Results)
}