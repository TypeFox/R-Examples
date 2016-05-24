ABautopacf <-
function(behavior,phaseX,v,lags){
  t1<-table(phaseX)
  tmax<-t1[names(t1)==v]
  start<-match(v,phaseX)
  end<-tmax+start-1
  tsx<-behavior[start:end]
  e=length(tsx)
  x=1:end
  x<-ts(x,start=1,end=e,deltat=1)
  graphics.off()
  pacf1<-pacf(tsx,lag.max=lags,plot=FALSE,na.action=na.pass,type=c("correlation"))
  b1<-pacf(tsx,lag.max=lags,plot=TRUE,na.action=na.pass)
  print(pacf1)
  print(b1)
  
}
