ABma <-
function(behavior,phaseX,v1){
 
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  t<-SMA(A,n=2)
  graphics.off()
  layout(rbind(1,2), heights=c(5,1))
  
  plot(A, lwd=2,type="l",col="blue", xlab="Time", ylab="behavior", bty='L',main="Moving Average  Plot" )
  
  lines(t,lwd=2,type="l",col="red",lty=2 )
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("MA","Behavior"), col = c("red", "blue"),lty=c(2,1), lwd = 1,ncol=2,bty ="n")
  behavior=c(A,NA)
  
 
  ma<-c(t,NA)
  endA=length(t)
  phase <- rep(v1, endA-1)
  phase <-c(phase,NA)
  ma<-ma[2:endA]
  ma<-c(ma,NA)
  transdat<-data.frame(ma,phase)
  a<-readline("Save results? (y/n) ")
  
  
  if (a=="y")
  {write.csv(transdat,file = tclvalue(tcl("tk_getSaveFile")),row.names=FALSE)}
  
  print(transdat)
}
