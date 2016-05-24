IQRline <-
function(behavior,phaseX,v){
  instruct<-"Click the mouse in the beginning of the phase you want the line in."
  print(instruct)
  abmedian<-tapply(behavior, phaseX,sd)
  meansd<-tapply(behavior, phaseX,mean)
  tphase<-table(phaseX)
  mlin<- tphase[names(tphase)==v]
  omedian<- abmedian[names(abmedian)==v]
  meansd1<-meansd[names(meansd)==v]
  
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v]
  startA<-match(v,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  q=quantile(A,na.rm=T)
  
  
  p75<-q[4]
  p25<-q[2]
  omedianu<-p75
  omedianb<-p25
  p2575<-c(omedianu,omedianb)
  print(p2575)
  l<-locator(1)
  mlin<-l$x+mlin+1
  segments(x0=l$x,x1=mlin,y0=omedianu,lwd=3,col="black")
  l<-locator(1)
  segments(x0=l$x,x1=mlin,y0=omedianb,lwd=3,col="black")
 u<-readline("accept line? (y/n) ")
  if (u=="n")
  {replayPlot(ab)}
  
    ab<-NULL
    ab<<-recordPlot()
}
