dohCluster <-function(X, peakList, refInd=0, 
    maxShift =100, acceptLostPeak=TRUE, verbose=TRUE){
  Y=X;
  peakListNew=peakList;
  
  if (verbose) startTime=proc.time();  
  refSpec=Y[refInd,];        
  for (tarInd in 1: nrow(X))
  if (tarInd!=refInd)
  {
      if (verbose) cat("\n aligning spectrum ",tarInd);
      targetSpec=Y[tarInd,];
      myPeakList=c(peakList[[refInd]],peakList[[tarInd]]);
      myPeakLabel=double(length(myPeakList));
      for (i in 1:length(peakList[[refInd]]) ) myPeakLabel[i]=1;
      startP=1;
      endP=length(targetSpec);
      res=hClustAlign(refSpec,targetSpec,myPeakList,myPeakLabel,startP,endP,
        maxShift=maxShift,acceptLostPeak=acceptLostPeak)
      Y[tarInd,]=res$tarSpec;

      peakListNew[[tarInd]]=
        res$PeakList[(length(peakList[[refInd]])+1):length(myPeakList)]
  }
  peakList=peakListNew;
  if (verbose){
    endTime=proc.time();
    cat("\n Alignment time: ",(endTime[3]-startTime[3])/60," minutes");
  }
  return(Y);
}
