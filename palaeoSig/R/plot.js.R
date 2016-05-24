plot.js<-function(x,names.v1,names.v2,...){
 EX<-x$EX
 v1<-x$v1
 v2<-x$v2
  plot(EX*v1,EX*v2, type="n", xlab="", ylab="", asp=1, xlim=range(c(EX*1.2,EX*-1.2) ), axes=F,...)
  polygon(EX*v1,EX*v2, col="grey90")
  points(EX*v1,EX*v2, pch=16,col=ifelse(round(EX,4)==round(max(EX),4),2,NA))
  lines((EX*v1)[round(EX,4)==round(max(EX),4)],(EX*v2)[round(EX,4)==round(max(EX),4)], pch=16,col=2, type="o")
  
  sig<-quantile(x$sim.ex, prob=.95)
  polygon(sig*v1,sig*v2, lty=2, border="blue")
  sapply(seq(0,max(EX), .1),function(r)  polygon(r*v1,r*v2, lty=3, border="grey60")) 
  arrows(0,0,sig*sqrt(.5)*c(-1,1,1,-1),sig*sqrt(.5)*c(-1,-1,1,1), length=.1)
  arrows(0,0,sig*c(-1,1,0,0),sig*c(0,0,-1,1), length=.1)

  text(sig*1.1*sqrt(.5)*c(-1,-1,1,1),sig*1.1*sqrt(.5)*c(-1,1,-1,1),labels=c(paste(names.v1[1],names.v2,sep="-"),paste(names.v1[2],names.v2,sep="-")))
  text(sig*c(-1.1,1.1,0,0),sig*c(0,0,-1.1,1.1), labels=c(names.v1,names.v2))
}
