ABautoacf <-
function(behavior,phaseX,v,l){
  t1<-table(phaseX)
  tmax<-t1[names(t1)==v]
  start<-match(v,phaseX)
  end<-tmax+start-1
  tsx<-behavior[start:end]
  e=length(tsx)
  x=1:end
  x<-ts(x,start=1,end=e,deltat=1)
  graphics.off()
  acf1<-acf(tsx,lag.max=l,plot=FALSE,na.action=na.pass,type=c("correlation"))
  b<-acf(tsx,lag.max=l,plot=TRUE,na.action=na.pass,type=c("correlation"))
  sig <- qnorm((1 + 0.95)/2)/sqrt(sum(!is.na(tsx)))
  b<-Box.test(tsx,lag=l, type="Ljung-Box")
   print(acf1)
  print(b)
  
 
}
