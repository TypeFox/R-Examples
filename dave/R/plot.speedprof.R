plot.speedprof <-
function(x,...) {
  o.speedprof<- x
  par(mfrow=c(1,1),omi=c(2,0,0,0))
# first comes matrixplot
  dimde<- o.speedprof$nrel
  maxord<- length(orders)
  mde<- o.speedprof$dmatrix
  timescale<- o.speedprof$timescale
  orders<- o.speedprof$orders
  capacity<- 0
  if(capacity == 0) capacity<- dimde
#  plot(c(-2,capacity),c(-2,capacity),type="n",xlab="",ylab="",ylim=c(capacity+1,-1),axes=FALSE,asp=1)
#  for (i in 1:dimde) for (j in 1:dimde){
#     points(i,j,cex=mde[i,j]/max(mde)*1.5,pch=1,col="black")
#     points(i,j,cex=mde[i,j]/max(mde)*1.5,pch=16,col="gray")
#  }
# first and last time step
#  text(rep(1,dimde),seq(1,dimde,1), timescale,pos=2,cex=0.6,offset=1.0)
#  text(seq(1,dimde,1),rep(1,dimde), timescale,pos=3,cex=0.6,srt=90,offset=1.2)
# now speedprofiles
  dd<- dimde-1
  dvec<- rep(1,dimde)
  jj<- 0
  linew<- seq(1,maxord,1)
  linew<- (linew*0.8)^1.5
# the exponent above causes spread of line width
  for (i in 1:dd) {
    jj<- jj+1
    j<- i+1
    dvec[jj]<- mde[i,j]
  }
  plot(c(min(timescale),max(timescale)),c(0,max(dvec)*1.1),type="n",xlab="Time scale",ylab="Rate of change",tcl=-0.3,mgp=c(2,0.5,0),cex.axis=0.8)
  colors<- rgb(0,0,0,seq(240,20,-(220/maxord+1)),maxColorValue=255)
  for (io in 1:maxord) {
    iom<- dimde-orders[io]
    dvec<- rep(1,iom)
    jj<- 0
    for (i in 1:iom) {
       ibeg<- 1+orders[io] 
       j<- i+orders[io]
       jj<- jj+1
       dvec[jj]<- mde[i,j]
    }
    dvec<- dvec/orders[io]
    lines(timescale[ibeg:dimde],dvec,lwd=linew[io],col=colors[io])
    lines(c(timescale[1], timescale[ibeg]),c(dvec[1],dvec[1]),lwd=linew[io],col=colors[io],lty="dotted")
 }
 legend("top",as.character(orders),title="Order (time steps involved)",lwd=linew,col=colors[1:maxord],bty="n",ncol=5,cex=0.8)

}
