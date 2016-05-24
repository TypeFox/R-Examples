SDAband <-
function(behavior,phaseX,v,bandX){
  
  abmedian<-tapply(behavior, phaseX,sd)
  meansd<-tapply(behavior, phaseX,mean)
  tphase<-table(phaseX)
  mlin<- tphase[names(tphase)==v]
  omedian<- abmedian[names(abmedian)==v]
  meansd1<-meansd[names(meansd)==v]
  omedianu<-omedian*bandX+meansd1
  omedianb<-meansd1-omedian*bandX
  
  l<-locator(1)
  mlin<-l$x+mlin+1
  segments(x0=l$x,x1=mlin,y0=omedianu,lwd=3,col="gray")
  l<-locator(1)
  segments(x0=l$x,x1=mlin,y0=omedianb,lwd=3,col="gray")
  
  u<-readline("accept line? (y/n) ")
  if (u=="n")
  {replayPlot(ab)}
  
    ab<-NULL
  ab<<-recordPlot()
  
}
