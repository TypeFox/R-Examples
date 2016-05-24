detectSpecPeaks <-function(X, nDivRange=128, scales=seq(1,16,2), 
    baselineThresh=50000, SNR.Th=-1, verbose=TRUE){
nFea=ncol(X);
nSamp=nrow(X);
noiseEsp=0.005;
if (SNR.Th<0) SNR.Th=max(scales) * 0.05; 
pList=NULL;
for (i in 1:nSamp){
  myPeakRes=NULL;
  mySpec=X[i,];
  for (k in 1:length(nDivRange)){
    divR=nDivRange[k];
    for (j in 1:(trunc(nFea/divR)-3)){
    startR=(j-1)*divR+1;
    if (startR>=nFea) startR=nFea;
    endR=(j+3)*divR;
    if (endR>nFea) endR=nFea;
    
    xRange=mySpec[startR:endR];
    xMean=mean(xRange);
    xMedian=median(xRange);
    if ((xMean==xMedian)||abs(xMean-xMedian)/((xMean+xMedian)*2)<noiseEsp){
      next;
    }
    else{
      peakInfo=
        peakDetectionCWT(mySpec[startR:endR],scales=scales,SNR.Th=SNR.Th)
      majorPeakInfo=peakInfo$majorPeakInfo
      if (length(majorPeakInfo$peakIndex)>0){
          myPeakRes=c(myPeakRes,majorPeakInfo$peakIndex+startR-1);
      }
    }
    }
  }
  pList[i]=list(myPeakRes);
  pList[[i]]=sort(unique(pList[[i]]));
  pList[[i]]=pList[[i]][which(X[i,pList[[i]]]>baselineThresh)]
  pList[[i]]=sort(pList[[i]]);
  if (verbose) cat("\n Spectrum ",i," has ",length(pList[[i]])," peaks");
}

return (pList)
}
