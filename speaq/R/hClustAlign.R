hClustAlign <-function(refSpec, tarSpec, peakList, peakLabel, startP, endP,
           distanceMethod="average", maxShift=0, acceptLostPeak=FALSE){
minpeakList=min(peakList)[1];
maxpeakList=max(peakList)[1];
startCheckP=startP+which.min(tarSpec[startP:(minpeakList-1)])[1]-1;
if (is.na(startCheckP)) startCheckP=startP;

endCheckP=maxpeakList+ which.min(tarSpec[(maxpeakList+1):endP])[1];
if (is.na(endCheckP)) endCheckP=endP;

if ((endCheckP-startCheckP)<2) {
  return (list(tarSpec=tarSpec,peakList=peakList));
}
adj=findShiftStepFFT(refSpec[startCheckP:endCheckP],
                          tarSpec[startCheckP:endCheckP],maxShift=maxShift);
if (adj$stepAdj!=0){
    if (acceptLostPeak) isImplementShift=TRUE 
    else isImplementShift=(adj$stepAdj<0&&adj$stepAdj+
      minpeakList >=startCheckP )||(adj$stepAdj>0&&adj$stepAdj+
      maxpeakList<=endCheckP);
    if (isImplementShift)
    {
      newTargetSpecRegion=doShift(tarSpec[startCheckP:endCheckP],adj$stepAdj);
      tarSpec[startCheckP:endCheckP]=newTargetSpecRegion;
      
      peakListTarget=which(peakLabel==0);
      peakList[peakListTarget]=peakList[peakListTarget]+adj$stepAdj;
      
        lostPeaks=which(peakList<=0);
        if (length(lostPeaks) >0){

          peakList=peakList[-lostPeaks];
          peakLabel=peakLabel[-lostPeaks];
      }
    }
}

if (length(peakList)<3) {return (list(tarSpec=tarSpec,peakList=peakList));}
hc=hclust(dist(peakList),method=distanceMethod)
clusterLabel=cutree(hc,h=hc$height[length(hc$height)-1]);
if (length(unique(clusterLabel))<2){ 
  return (list(tarSpec=tarSpec,peakList=peakList));
  }

labelID1=which(clusterLabel==1);
subData1=peakList[labelID1];
subLabel1=peakLabel[labelID1];

labelID2=which(clusterLabel==2);
subData2=peakList[labelID2];
subLabel2=peakLabel[labelID2];
maxsubData1=max(subData1)[1];
minsubData2=min(subData2)[1];

if (maxsubData1<minsubData2){
    endP1=maxsubData1+which.min(tarSpec[(maxsubData1+1) :(minsubData2-1)])[1];
  if (is.na(endP1)) endP1=maxsubData1;
    startP2=endP1+1;
    if (length(unique(subLabel1))>1){
        res=hClustAlign(refSpec,tarSpec,subData1,subLabel1,startP,endP1,
          distanceMethod=distanceMethod,maxShift=maxShift,
          acceptLostPeak=acceptLostPeak);
        tarSpec=res$tarSpec;
        peakList[labelID1]=res$peakList;
    }        
    if (length(unique(subLabel2))>1){
        res=hClustAlign(refSpec,tarSpec,subData2,subLabel2,startP2,endP,
          distanceMethod=distanceMethod,maxShift=maxShift,
          acceptLostPeak=acceptLostPeak);
        tarSpec=res$tarSpec;
        peakList[labelID2]=res$peakList;
    }
}else{        
    maxsubData2=max(subData2)[1];
    minsubData1=min(subData1)[1];
    endP2=maxsubData2+which.min(tarSpec[(maxsubData2+1) :(minsubData1-1)])[1];
   if (is.na(endP2)) endP2=maxsubData2;
    startP1=endP2+1;
    if (length(unique(subLabel2))>1){
        res=hClustAlign(refSpec,tarSpec,subData2,subLabel2,startP,endP2,
          distanceMethod=distanceMethod,maxShift=maxShift,
          acceptLostPeak=acceptLostPeak);
        tarSpec=res$tarSpec;
        peakList[labelID2]=res$peakList;
    }    
    if (length(unique(subLabel1))>1){
        res=hClustAlign(refSpec,tarSpec,subData1,subLabel1,startP1,endP,
          distanceMethod=distanceMethod,maxShift=maxShift,
          acceptLostPeak=acceptLostPeak);
        tarSpec=res$tarSpec;
        peakList[labelID1]=res$peakList;
    }        
}    
return (list(tarSpec=tarSpec,peakList=peakList));
}