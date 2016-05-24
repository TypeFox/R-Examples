drawBW <-function(BW, perc, X, startP=-1, endP=-1, groupLabel=NULL,
     highBound=-1, lowBound=-1, nAxisPos=4, offside=0){
  if (startP==-1) startP=1;
  if (endP==-1) endP=ncol(X);
  GraphRange<-c(startP:endP); 
  op=par(mfrow=c(2,1))
      
  plot(BW[GraphRange],ylim=c(min(c(BW[GraphRange]),perc[GraphRange]),
    max(c(BW[GraphRange]),perc[GraphRange])),type="n", ylab="BW",
  xlab="index", xaxt="n")
  tempVal =trunc(length(GraphRange)/nAxisPos);
  xPos=c(0:nAxisPos) * tempVal; 
  axis(1,at=xPos,labels=xPos+startP+offside);
  lines(BW[GraphRange],col= "blue")
  lines(perc[GraphRange],col= "black")
  drawSpec(X,xlab="index",ylab="intensity",startP=startP,endP=endP,groupLabel=groupLabel,
    offside=offside,highBound=highBound,lowBound=lowBound)
}
