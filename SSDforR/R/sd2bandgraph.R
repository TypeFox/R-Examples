sd2bandgraph <-
function(behavior,phaseX,v1,ABxlab,ABylab,ABmain){
  t1<-table(phaseX)
  tmaxA<-t1[names(t1)==v1]
  startA<-match(v1,phaseX)
  endA<-tmaxA+startA-1
  A<-behavior[startA:endA]
  
  meanA=mean(A,na.rm=T)
  sdA=sd(A,na.rm=T)
  x1=(c(seq(1:tmaxA)))
  SDabove<-meanA+sdA*2
  SDbelow<-meanA-sdA*2
  
  maxy=which.max(behavior)
  max<-behavior[maxy]
  #numx<-sum(!is.na(behavior))+1
  numx=length(x1)+1
  max=SDabove+2
  min=SDbelow-2
  
  graphics.off()
  layout(rbind(1,2), heights=c(6,1))
  plot(A,ylim=c(min,max),xlim=c(0,numx),type="o",xlab=ABxlab, ylab=ABylab, main=ABmain,bty='L')
  abline(h=SDabove,col="blue")
  abline(h=SDbelow,col="red")
  abline(h=meanA,col="orange")
  par(mar=c(1, 1, 1, 1))
  plot.new()
  legend("center", c("+2sd","mean","-2sd"), col = c("blue","orange", "red"), lwd = 1,ncol=3,bty ="n")
  sd1<-c("SD=",round(sdA,2))
  psdu<-c("+2sd=",round(SDabove,2))
  pmean<-c("mean=",round(meanA,2))
  psdb<-c("-2SD=",round(SDbelow,2))
  print(sd1)
  print(psdu)
  print(pmean)
  print(psdb)
  ab<-NULL
  ab<<-recordPlot()
}
