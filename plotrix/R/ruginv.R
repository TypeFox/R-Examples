ruginv<-function(x,ticksize=0.03,side=1,lwd=0.5,col=par("fg"),
 col.ticks="white",quiet=getOption("warn")<0,...) {

 x<-as.vector(x)
 ok=is.finite(x)
 x<-x[ok]
 if(!quiet) {
  u<-par("usr")
  if(side%%2 == 1) {
   if(par("xlog")) u<-10^u[1:2]
   else u<-u[1:2]
  }
  else {
   if(par("ylog")) u<-10^u[3:4]
   else u<-u[3:4]
  }
  if(any(x < u[1] | x > u[2])) warning("Some values will be clipped")
 }
 u<-par("usr")
 par("pin")
 if(ticksize < 0.5) tic<-min(diff(u[3:4]),diff(u[1:2]))*ticksize
 else tic<-ifelse(side%%2 == 1,diff(u[3:4]),diff(u[1:2]))*ticksize
 if(ticksize < 0) opar<-par(xpd=TRUE)
 switch(as.character(side),
  "1"=polygon(u[c(1,2,2,1,1)],u[3]+c(0,0,tic,tic,0),col=col,border=NA,...),
  "2"=polygon(u[1]+c(0,0,tic,tic,0),u[c(3,4,4,3,3)],col=col,border=NA,...),
  "3"=polygon(u[c(1,2,2,1,1)],u[4]+c(0,0,-tic,-tic,0),col=col,border=NA,...),
  "4"=polygon(u[2]+c(0,0,-tic,-tic,0),u[c(3,4,4,3,3)],col=col,border=NA,...))
 switch(as.character(side),
  "1"=sapply(x,function(z) lines(c(z,z),u[3]+c(0,tic),col=col.ticks,lwd=lwd)),
  "2"=sapply(x,function(z) lines(u[1]+c(0,tic),c(z,z),col=col.ticks,lwd=lwd)),
  "3"=sapply(x,function(z) lines(c(z,z),u[4]+c(0,-tic),col=col.ticks,lwd=lwd)),
  "4"=sapply(x,function(z) lines(u[2]+c(0,-tic),c(z,z),col=col.ticks,lwd=lwd)))
 if(ticksize < 0) par(opar)
 invisible(x)
}

