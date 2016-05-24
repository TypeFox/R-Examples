IQRbandgraph <-
function(behavior,phaseX,v1,ABxlab,ABylab,ABmain){
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  q=quantile(A,na.rm=T)
  
  x1=(c(seq(1:tmaxA)))
  p75<-q[4]
  p25<-q[2]
  medianA<-q[3]
  maxy=which.max(behavior)
  max<-behavior[maxy]
  #numx<-sum(!is.na(behavior))+1
  numx=length(x1)+1
  graphics.off()
  layout(rbind(1,2), heights=c(6,1))
  plot(A,ylim=c(0,max),xlim=c(0,numx),type="o",xlab=ABxlab, ylab=ABylab, main=ABmain,bty='L')
  abline(h=p75,col="blue")
  abline(h=p25,col="red")
  abline(h=medianA,col="orange")
  
  
  psdu<-c(round(p75,2))
  pmean<-c(round(medianA,2))
  psdb<-c(round(p25,2))
  iqr=q[4]-q[2]
  iqrp=c("IQR=",round(iqr,2))
  
  tprint=c(psdu,pmean,psdb)
  print(tprint)
  ab<-NULL
  ab<<-recordPlot()
}
