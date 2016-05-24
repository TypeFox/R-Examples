Arimama <-
function(behavior,phaseX,v,m){
  
t1<-table(phaseX)
tmax<-t1[names(t1)==v]
start<-match(v,phaseX)
end<-tmax+start-1
tsx<-behavior[start:end]
e=length(tsx)
x=1:end
x<-ts(x,start=1,end=e,deltat=1)
tsma<-SMA(tsx,n=m)
graphics.off()
par(mfrow=c(1,2))
plot.ts(tsx,xlab="time", ylab="behavior")
plot.ts(tsma,xlab="time", ylab="Moving Average of behavior")
}
