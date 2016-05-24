plot.mxplot <-
function(x,...,capacity=100) {
  o.mxplot<- x
  par(mai=c(1,1,1,1))
  par(mar=c(4,1,1,1))
  mc<- o.mxplot$order
  com<- o.mxplot$mmatrix
  lev<- o.mxplot$levels
  com<- com-min(com)
  com<- com/max(com)
  plot(c(0.5,capacity),c(0.5,capacity),type="n",xlab="",ylab="",ylim=c(capacity,0.5),axes=FALSE,asp=1)
  i<- rep(seq(1,mc,1),mc)
  j<- rep(1:mc,rep(mc,mc))
  pscal<- (100/capacity)*0.9
  points(i,j,pch=16,cex=com*pscal,col="gray")
  points(i,j,pch=1,cex=com*pscal)
  i<- seq(1,mc,1)
  text(rep(1,mc),i,lev,pos=2,cex=0.4,offset=0.6)
  text(i,c(rep(mc,mc)),lev,pos=1,cex=0.4,offset=1.5,srt=90)
  }
