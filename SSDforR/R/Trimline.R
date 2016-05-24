Trimline <-
function(behavior,phaseX,v){
  writeLines("-------------------------------------------------------------------------------------")
  writeLines("Click the mouse in the beginning of the phase you want the line in")
  writeLines("-------------------------------------------------------------------------------------")
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
  trimmean=mean(A,trim=.1,na.rm=T)
  trm<-c("10% Trim Mean=",round(trimmean,3))
  print(trm)
  l<-locator(1)
  mlin<-l$x+mlin
  segments(x0=l$x,x1=mlin,y0=trimmean,col="red",lty=2,lwd=3)
  
  u<-readline("accept line? (y/n) ")
  if (u=="n")
  {replayPlot(ab)}
  
 ab<-NULL 
  ab<<-recordPlot()
  
 
}
