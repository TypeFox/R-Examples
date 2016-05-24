"waitGUI"<-function(i,itr){


temp<-function(z){
plot(1,1,ylim=c(0,1),col="white",xlim=c(0,100),ylab=NA,xlab=NA,xaxt="n",yaxt="n",axes=FALSE)
segments(0,0.47,100,0.47,lwd=5)
segments(0,0.53,100,0.53,lwd=5)
segments(0,0.47,0,0.53,lwd=5)
segments(100,0.47,100,0.53,lwd=5)
a<-length(c(0:z))
if(a >= 101) a<-100
segments(c(0:a),rep(0.47,a),c(0:a),rep(0.53,a),lwd=rep(5,a))
text(2,0.43,"0 %")
text(99,0.43,"100 %")
text(50,0.57,"Please, wait while working")}
if (i==1) z<<-0
if (i==1) temp(z)
if (itr<=99) z<<-c((i*100)/itr)
if (itr<=99) temp(z)
if (itr>=100) z<<- round(c((i*100)/itr),digits=0)
if (itr>=100) temp(z)
if (z==100) z<<-99
if (z==100)   temp(z)

}
