CalculateHUM_Ex<-function(data,indexF,indexClass,allLabel,amountL)
{
  #library(Rcpp)
  #library(gtools)
    
  dataEach=NULL
 
  class.index=NULL
  for(i in 1:length(allLabel))
  {
      vrem=which(data[,indexClass]==allLabel[i])
      dataEach=c(dataEach,list(data[vrem,indexF,drop = FALSE]))
      class.index=c(class.index,list(vrem))
  }
  
  
  indexLabel<-combn(allLabel,amountL)
  output<-matrix(ncol=(length(indexF)+amountL),nrow=ncol(indexLabel))
 
  seqAll=NULL
  
  #cycle for different label combinations
  for(j in 1:ncol(indexLabel))
  {
    indexUnion=NULL
    indexL=NULL
    
    for(i in 1:amountL)
    {
      v.class=which(allLabel==indexLabel[i,j])
      indexL=c(indexL,v.class)
      output[j,i]<-indexLabel[i,j]
      indexUnion=union(indexUnion,class.index[[v.class]])
    }
  
   seqMax=NULL
   seq=gtools::permutations(amountL,amountL,1:amountL)
    
  for(i in 1:length(indexF))
  {
    
    s_data=NULL
    dataV=data[,indexF[i]]
     
    prodValue=1
    for (k in 1:amountL)
    {
      vrem=sort(dataEach[[indexL[k]]][,i])
    
      s_data=c(s_data,list(vrem))
      prodValue = prodValue*length(vrem)
    }
    
    #claculate the threshold values for plot of 2D ROC and 3D ROC
    thresholds <- sort(unique(dataV[indexUnion]))
    thresholds=(c(-Inf, thresholds) + c(thresholds, +Inf))/2
    
    out=CalcGene(s_data,seq, prodValue,thresholds)
    
    output[j,(amountL+i)]<-out$HUM
    seqMax=cbind(seqMax,out$seq)
    
  }
  colnames(seqMax)=names(data[,indexF,drop = FALSE])  
  
  seqAll=c(seqAll,list(seqMax))
  } 
  
    name<-NULL
    for(i in 1:amountL)
    {
      name<-c(name,paste("Diagnosis",i,sep=""))
    }
    colnames(output)<-c(name,indexF)
  

  return(list(HUM=output,seq=seqAll))
}
