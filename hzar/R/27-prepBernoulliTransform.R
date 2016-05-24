

## mkn.charLocTwinV
tVec.make <- function(loc,locDist,char,charLoc){
  keepVals=!is.na(char);
  char=char[keepVals]
  charLoc=charLoc[keepVals];
##   print(loc);
##   print(charLoc[1:10]);
##   print(char[1:10]);
  locF <- factor(loc,sort(unique(c(as.character(loc),
                                   as.character(charLoc)))))
   ##  mkn.fixFactor(loc,charLoc);         
## print(locF);
  charLocF <- factor(charLoc,levels(locF));
  distOfLoc <- locDist;names(distOfLoc) <- locF;
  distOfLoc[sapply(levels(locF), function(x) which(x==locF)[[1]])]->distOfLoc;
  return(list(dLoc=distOfLoc,
              lChar=data.frame(loc=charLocF,chData=char)));
}

mkn.fixFactor <- function(x,levelsOfX=levels(x)){
  return(factor(x,sort(unique(c(as.character(x),
                                as.character(levelsOfX))))));
}

bLLdist <- function(series,totals){
  hS <- ifelse(series>0,series*log(totals/series),0);

  return(sum(hS));
}

wMax<-function(x)which(x==max(x))



bLLSumLocDistFreq <- function(tVec,cutValues,
                    values=tVec$lChar$chData,
                    locations=tVec$lChar$loc){
  cV=sort(unique(cutValues));
  pL=lapply(cV,function(x) locations[values<x]);
  pL=c(pL,list(locations));
  
  cL=lapply(pL,summary);
  sL=c(cL[1],lapply(2:length(cL),function(x)cL[[x]]-cL[[x-1]]))
  mL=matrix(unlist(sL),ncol=length(cL),byrow=FALSE)
  mP=mL/apply(mL,1,sum)
  
  return(sum(ifelse(mP>0,
                    mP*log(rep(apply(mP,2,sum),
                               each=nrow(mP)
                               ) / mP
                           ),
                    0)) );
  ## bLLdist(mP,rep(apply(mP,2,sum),each=nrow(mP)));
}

midpoints<-function(x) {
    y<-as.numeric(sort(unique(x)));
    return((y[y>min(y)]+y[y<max(y)])/2);
  }


bLLSumDistSF <- function(tVec,
                    values=tVec$lChar$chData,
                    locations=tVec$lChar$loc){
  mVal=midpoints(values);
  data.frame(char=mVal,info=sapply(mVal,bLLSumLocDistFreq,tVec=tVec,values=values,locations=locations))
}
bLLSLDSFchoice <- function(tVec,
                    values=tVec$lChar$chData,
                    locations=tVec$lChar$loc){
  mVal=midpoints(values);
  return(mVal[which.min(sapply(mVal,bLLSumLocDistFreq,tVec=tVec,values=values,locations=locations))]);
}


pullMinScore <- function(bFrame) bFrame$char[wMax(-bFrame$info)]


hzar.dBernoulli.LL <- function(values,locations,getMax=FALSE,getProbs=FALSE){
  if(getProbs){
    cV <- bLLSLDSFchoice(NULL,values=values,locations=locations);
  pL<-lapply(unique(locations),function(x) as.numeric(values[locations==x]));
    nSamples<-as.numeric(lapply(pL,length));
    nLow<-as.numeric(lapply(pL,function(x,cutValue) sum(x<cutValue), cV) )
    return(data.frame(locID=unique(locations),
                      pObs=nLow/nSamples,
                      nLow=nLow,
                      nEff=nSamples,
                      cutValue=cV));
    
    
  }
  if(getMax){
    return(bLLSLDSFchoice(NULL,values=values,locations=locations))
  }
  return(bLLSumDistSF(NULL,values=values,locations=locations)$info);
} 
tVec.dBernoulli.LL <- function(tVec,
                               values=tVec$lChar$chData,
                               locations=tVec$lChar$loc,
                               getMax=FALSE,getProbs=FALSE){
  if(getProbs){
    cV <- bLLSLDSFchoice(tVec);
    pL<-lapply(unique(locations),function(x) as.numeric(values[locations==x]));
    nSamples<-as.numeric(lapply(pL,length));
    nLow<-as.numeric(lapply(pL,function(x,cutValue) sum(x<cutValue), cV) )
    return(data.frame(locID=unique(locations),
                      pObs=nLow/nSamples,
                      nLow=nLow,
                      nEff=nSamples,
                      cutValue=cV));
    
    
  }
  if(getMax){
    return(bLLSLDSFchoice(tVec));
  }
  return(bLLSumDistSF(tVec)$info);
} 
old.dBernoulli.LL <- function(values,locations,getMax=FALSE,getProbs=FALSE){
  
  mVal<-midpoints(values);
  bernoulli.LL <- function(nSuccess,nTotal){
    nFailures<-nTotal-nSuccess;
    hSuccess<-ifelse(nSuccess>0,nSuccess*log(nTotal/nSuccess),0);
    hFailure<-ifelse(nFailures>0,nFailures*log(nTotal/nFailures),0);
    return(hSuccess+hFailure);
  }
  pL<-lapply(unique(locations),function(x) as.numeric(values[locations==x]));
  bLLSum<- function(junk,procList){
    nSamples=as.numeric(lapply(procList,length))
    return(sum(bernoulli.LL(as.numeric(lapply(procList,
                                              function(x,cutValue)
                                              sum(x<cutValue),
                                              cutValue=junk)),
                            nSamples)));
  }
  if(getProbs){
    junk<-as.numeric(lapply(mVal, bLLSum, pL ));
    cutValue<-mVal[junk==max(junk)][1];
    nSamples<-as.numeric(lapply(pL,length));
    nLow<-as.numeric(lapply(pL,function(x,cutValue) sum(x<cutValue), cutValue) )
    return(data.frame(locID=unique(locations),
                      pObs=nLow/nSamples,
                      nLow=nLow,
                      nEff=nSamples,
                      cutValue=cutValue));
  }
  if(getMax){
    junk<-as.numeric(lapply(mVal, bLLSum, pL ));
    return(mVal[junk==max(junk)][1]);
  }
  return(as.numeric(lapply(mVal, bLLSum, pL )));
}

hzar.makeTraitObsData <- function(distOfLocation,locationOfValue,values){
  probs <- hzar.dBernoulli.LL(values,locationOfValue,getProbs=TRUE)
  return(hzar.doMolecularData1DPops(distance=distOfLocation[probs$locID],
                                    pObs=probs$pObs,
                                    nEff=probs$nEff,
                                    siteID=probs$locID));
}

tVec.makeTraitObsData <- function(tVec){
  probs <- tVec.dBernoulli.LL(tVec,getProbs=TRUE)
  return(hzar.doMolecularData1DPops(distance=tVec$dLoc[probs$locID],
                                    pObs=probs$pObs,
                                    nEff=probs$nEff,
                                    siteID=probs$locID));
}


hzar.doMorphoSets <- function(traitNames, tDist, tDLocCol, tDDistCol, tValues, tVLocCol){

  res <- lapply(traitNames,
                function(id) {
                  return(tVec.makeTraitObsData(tVec.make(loc=tDist[[tDLocCol]],
                                                                locDist=tDist[[tDDistCol]],
                                                                char=tValues[[id]],
                                                                charLoc=tValues[[tVLocCol]])));
                  }
                )
  names(res) <- traitNames;
  return(res);
}

old.doMorphoSets <- function(traitNames, tDist, tDLocCol, tDDistCol, tValues, tVLocCol){
  distOfLocation <- tDist[[tDDistCol]];
  names(distOfLocation) <- tDist[[tDLocCol]];
  distOfLocation[sapply(levels(tDist[[tDLocCol]]),
                        function(x,ids) which(ids==x)[[1]],
                        tDist[[tDLocCol]])]->distOfLocation;

  res <- lapply(traitNames,
                function(id) {
                  tVal <- tValues[!is.na(tValues[[id]]),c(tVLocCol,id)]
                  
                  return(hzar.makeTraitObsData(distOfLocation,
                                        locationOfValue=tVal[[tVLocCol]],
                                        values=tVal[[id]]));
                  }
                )
  names(res) <- traitNames;
  return(res);
}
