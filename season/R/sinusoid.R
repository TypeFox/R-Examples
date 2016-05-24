sinusoid<-function(amplitude,frequency,phase,...){
  time<-seq(0,2*pi,pi/1000);
  sinusoid<-amplitude*cos(time*frequency-phase)
  par(las=1)
  plot(time,sinusoid,type='l',xaxt='n',xlab='Time (radians)',ylab='',...)
  lines(range(time),c(0,0),lty=2)
  axis(side=1,at=c(0,pi,2*pi),font=5,labels=c('0','p','2p'))
}
#sinusoid(amplitude=1,frequency=1,phase=1)