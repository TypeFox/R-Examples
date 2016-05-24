findSegPeakList <-function(peakList, startP, endP){
  res=0;  
  for (i in 1:length(peakList)){    
    if (peakList[i]>startP&&peakList[i]<endP){    
      res=c(res,peakList[i]-startP+1);
    }    
  }
  return(res[-1])
}
