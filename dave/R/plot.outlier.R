plot.outlier <-
function(x,...) {
     o.outlier<- x
     neighb<- o.outlier$neigh.dist
     dmm<- o.outlier$olddim[1]
     par(omi=c(1,0,1.5,0))
     hist(neighb,xlab="Nearest neighbour distance",main="",col=gray(0.5),labels=T,ylim=c(0,dmm*0.6))
# ordination
     x.c<- o.outlier$pco.points[,1]
     y.c<- o.outlier$pco.points[,2]
     thresh<- o.outlier$threshold
     par(omi=c(0,0,0,0))
     plot(x.c,y.c,asp=1,type="n",cex.axis=0.7,cex.lab=1.0,tcl=-0.3,mgp=c(2,0.5,0))
     points(x.c[neighb <= thresh],y.c[neighb <= thresh],pch=16,cex=0.6,col=gray(0.8))
     points(x.c[neighb > thresh],y.c[neighb > thresh],pch=16,cex=0.6)
     abline(h=0,v=0,lwd=1.0,col="gray")
     legend("bottomleft",c("outliers","other data points"),pch=c(16,16),col=c(gray(0.0),gray(0.8)),bty="n",cex=0.8)
#     legend("topleft","A",cex=1.5,bty="n",inset=c(-0.05,-0.02))
   }
