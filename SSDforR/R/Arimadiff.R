Arimadiff <-
function(behavior,phaseX,v,d){
  t1<-table(phaseX)
  tmax<-t1[names(t1)==v]
  start<-match(v,phaseX)
  end<-tmax+start-1
  tsx<-behavior[start:end]
  e=length(tsx)
  x=1:end
  x<-ts(x,start=1,end=e,deltat=1)
  tsdiff<-diff(tsx,differences=d)
  graphics.off()
  par(mfrow=c(1,2))
  plot.ts(tsx,xlab="time", ylab="behavior")
  plot.ts(tsdiff,xlab="time", ylab="difference of behavior")
}
