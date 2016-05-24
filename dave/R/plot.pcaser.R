plot.pcaser <-
function(x,...,lines=TRUE) {
  o.pcaser<- x
  par(mfrow=c(1,1),omi=c(0.5,0,0,0),mgp=c(1.5,0.5,0),pty="s")
#
  plot(o.pcaser$scores[,1],o.pcaser$scores[,2],type="n",xlab="PCA axis 1",ylab="PCA axis 2",asp=1,cex.axis=0.8,cex.lab=0.8,tcl=-0.2)
  abline(h=0,v=0,lwd=1.0,col="gray")
  ser<- o.pcaser$plotlab
  points(o.pcaser$scores[,1],o.pcaser$scores[,2],pch=ser,cex=1.0)
  symbols<- as.integer(levels(as.factor(ser)))
  legend("topleft",levels(o.pcaser$plotlabels),pch=symbols,cex=0.8,bty="n",ncol=3)
  if(lines == TRUE) {
      for (i in 1:o.pcaser$nser) {
          for (j in 1:o.pcaser$nrel) {
             lines(o.pcaser$scores[ser==i,1],o.pcaser$scores[ser==i,2],lwd=0.5,col="gray")
          }
      }
   }
 }
