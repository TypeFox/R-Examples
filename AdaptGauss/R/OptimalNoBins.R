OptimalNoBins <-
function(Data){
# function [OptimalNrOfBins] = OptNrOfBins(data); 
# % [OptimalNrOfBins] = OptNrOfBins(data)
# %
# % DESCRIPTION
# % Berechung der optimalen Anzahl von Bins fuer ein Histogramm
# % nach Keating/Scott 99
# % INPUT
# % data               die Daten
# % OUTPUT
# % OptimalNrOfBins   die bestmoegliche ANzahl von Bins, minimal jedoch 10
# %                   Verwendung fuer hist(data,OptimalNrOfBins);
  #Anzahl vorhandene Daten
    if(is.matrix(Data)) nData <- colSums(!is.nan(Data))
    if(is.vector(Data)) nData <- sum(!is.nan(Data))
    
    prctile_hlp=function(x,p){
#   matlab:
#   Y = prctile(X,p) returns percentiles of the values in X. 
#   p is a scalar or a vector of percent values. When X is a 
#   vector, Y is the same size as p and Y(i) contains the p(i)th 
#   percentile. When X is a matrix, the ith row of Y contains the 
#   p(i)th percentiles of each column of X. For N-dimensional arrays,
#   prctile operates along the first nonsingleton dimension of X.  
  if(length(p)==1){  
            if(p>1){p=p/100}
            
  }

  if(is.matrix(x) && ncol(x)>1){
    cols<-ncol(x)
    quants<-matrix(0,nrow=length(p),ncol=cols)
    for(i in 1:cols){
      quants[,i]<-quantile(x[,i],probs=p,type=5,na.rm=TRUE)
    }
  }else{
    quants<-quantile(x,p,type=5,na.rm=TRUE)
  }
  return(quants)
} 
    
    
    if(nData<1){
      optNrOfBins<-0
    }else{    
      sigma<-sd(Data,na.rm=TRUE)    
      p<-prctile_hlp(Data,c(0.25,0.75))
      interquartilRange<-p[2]-p[1]     
      sigmaSir<-min(sigma,interquartilRange/1.349)
      optBinWidth<-3.49*sigmaSir/(nData)^(1/3)
      if(optBinWidth>0){
        optNrOfBins<-max(ceiling((max(Data,na.rm=TRUE)-min(Data,na.rm=TRUE))/optBinWidth),10)
      }else{
        optNrOfBins<-10
      }
    }                  
 return (optNrOfBins) 
 }