dohClusterCustommedSegments <-function(X, peakList, refInd, maxShift,
     acceptLostPeak=TRUE, segmentInfoMat, minSegSize=128,verbose=TRUE){
  if (!is.matrix(segmentInfoMat)) {
    cat("ERROR! segmentInfoMat must be in a matrix format.")
    return
  }

  if (segmentInfoMat[1,1]>minSegSize) 
    mysegments=c(1,segmentInfoMat[1,1]-1,0,0,0) else mysegments=c();
  i = 0;
  if (nrow(segmentInfoMat)>1){   
   for (i in 1:(nrow(segmentInfoMat)-1)){
    mysegments=c(mysegments,c(segmentInfoMat[i,]))
    if (segmentInfoMat[i+1,1]-segmentInfoMat[i,2]>minSegSize)
      mysegments=c(mysegments,c(segmentInfoMat[i,2]+1,segmentInfoMat[i+1,1]-1,0,0,0))
   }
  }
  mysegments=c(mysegments,c(segmentInfoMat[i+1,]))  

  if (ncol(X)-segmentInfoMat[i+1,2]-1>minSegSize) 
    mysegments=c(mysegments,c(segmentInfoMat[i+1,2]+1,ncol(X),0,0,0))

  mysegments=matrix(data=mysegments,nrow=length(mysegments)/5,
              ncol=5,byrow=TRUE,dimnames=NULL)
  
  mysegments[which(mysegments[,4]==0),4]=refInd;
  mysegments[which(mysegments[,5]==0),5]=maxShift;
  
  if (sum(mysegments[,3]!=0)==0){
    cat("\n No segments are set for alignment! Please set 
        at least 1 values in columnn 3 in segmentInfoMat matrix be 1 ")
    return;
  }

  Y=X;

  for (i in 1:nrow(mysegments))
  if (mysegments[i,3]!=0)
  {
    if (verbose)
    cat("\n Doing alignment a segment from ",
        mysegments[i,1]," to ",mysegments[i,2]," ...");

    segmentpeakList=peakList;
    for (j in 1:length(peakList)){
      segmentpeakList[[j]]=
          findSegPeakList(peakList[[j]],mysegments[i,1],mysegments[i,2]);
    }    
    Y[,c(mysegments[i,1]:mysegments[i,2])]=
    dohCluster(X[,c(mysegments[i,1]:mysegments[i,2])],peakList=segmentpeakList,
        refInd=mysegments[i,4],maxShift =mysegments[i,5],
        acceptLostPeak=acceptLostPeak, verbose=verbose);    
  }else{
    if (verbose)
      cat("\n The segment ",
      mysegments[i,1],"-",mysegments[i,2], " is not aligned");
  }  
 return(Y)
}