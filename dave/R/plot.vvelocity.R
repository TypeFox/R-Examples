plot.vvelocity <-
function(x,...,tlabs=c(1),scal=1) {
  o.vvelocity<- x
  pveg<- o.vvelocity$pveg
  timescale<- o.vvelocity$timescale
  y<- o.vvelocity$y
# speed graph
# -----------
  ntim <- length(pveg[,1])                                   # number of time steps
  time<- as.double(timescale)                                # the time vector
  y<- 0.5
  dmat<- dist(pveg^y,method="euclidean",diag=T,upper=T)      # distance matrix of time steps
  dmat<- as.matrix(dmat)
  size<- rep(0,ntim)
  for (i in 2:ntim) size[i]<- dmat[i-1,i]                    # taking elements below the diagonal of dmat
  size[1]<- size[2]                                          # first element the same as second
  out<- pco(dmat)
  par(lwd=0.5,mfrow=c(1,1))
  plot(out$points[,1]*scal,out$points[,2]*scal,asp=1,cex.axis=0.8,cex.lab=0.8,tcl=-0.2,ylab="PCO axis 2",xlab="PCO axis 1",type="n",frame.plot=T,tcl=-0.3,mgp=c(1.5,0.5,0))
  points(out$points[,1],out$points[,2],pch=21,bg="gray",cex=size*3)
  text(out$points[tlabs,1],out$points[tlabs,2],time[tlabs],pos=c(4,3,4,3,2,2,2),offset=1.5,cex=0.8)
  abline(h=0,v=0,col="gray")
# legend("topleft",expression(paste("PCoA, x' = ",x^y)),bty="n",cex=0.8)
# acceleration graph
# ------------------
  for (i in 2:ntim) size[i]<- dmat[i-1,i] 
  accel<- rep(0,ntim)
  shade<- rep("gray",ntim)
  for (i in 2:ntim){ 
    accel[i]<- size[i]-size[i-1]
    if (accel[i] <= 0.0) shade[i]<- "white"
  }
  par(lwd=0.5,mfrow=c(1,1))
  plot(out$points[,1]*scal,out$points[,2]*scal,asp=1,cex.axis=0.8,cex.lab=0.8,tcl=-0.2,ylab="PCO axis 2",xlab="PCO axis 1",type="n",frame.plot=T,tcl=-0.3,mgp=c(1.5,0.5,0))
  points(out$points[,1],out$points[,2],pch=21,bg=shade,cex=abs(accel)*20)
#  text(out$points[c(1,15,48,60,100,122,ntim),1],out$points[c(1,15,48,60,100,122,ntim),2],
#+ time[c(1,15,48,60,100,122,ntim)],pos=c(4,3,4,3,2,2,2),offset=1.5,cex=0.8)
  text(out$points[tlabs,1],out$points[tlabs,2],time[tlabs],pos=c(4,3,4,3,2,2,2),offset=1.5,cex=0.8)
  abline(h=0,v=0,col="gray")
  legend("bottomleft",c("Acceleration, positive","Acceleration, negative"),pch=21,pt.bg=c("gray","white"),pt.cex=1.5,bty="n",cex=0.8)
# scale speed profiles
# --------------------
  par(mfrow=c(4,1),mar=c(3,4,1,4))
# first graph with x'=x^1.0
  dmat<- dist(pveg^1.0,method="euclidean",diag=T,upper=T)    # distance matrix of time steps
  dmat<- as.matrix(dmat)
  size<- rep(0,ntim)
  for (i in 2:ntim) size[i]<- dmat[i-1,i]                    # taking elements below the diagonal of dmat
  size[1]<- size[2]   
  plot(time[2:ntim],size[2:ntim],xlab="",ylab="Rate of change",type="n",cex.axis=1.0,cex.lab=1.0,tcl=-0.2,mgp=c(2.0,0.5,0))
  lines(time[2:ntim],size[2:ntim],lwd=2.0,col=gray(0.6))                  
  legend("top",expression(paste(italic("x' = "),italic(x)^1.0)),bty="n",cex=1.2)
legend("topright","(a)",cex=1.6,bty="n")
# second graph with x'=x^0.25
  dmat<- dist(pveg^0.25,method="euclidean",diag=T,upper=T)    # distance matrix of time steps
  dmat<- as.matrix(dmat)
  size<- rep(0,ntim)
  for (i in 2:ntim) size[i]<- dmat[i-1,i]                    # taking elements below the diagonal of dmat
  size[1]<- size[2]   
  plot(time[2:ntim],size[2:ntim],xlab=" ",ylab="Rate of change",type="n",cex.axis=1.0,cex.lab=1.0,tcl=-0.3,mgp=c(2.0,0.5,0))
  lines(time[2:ntim],size[2:ntim],lwd=1.2,col=gray(0.6))                  
  legend("top",expression(paste(italic("x' = "),italic(x)^0.25)),bty="n",cex=1.2)
legend("topright","(b)",cex=1.6,bty="n")
# third graph with x'=x^0.05
  dmat<- dist(pveg^0.05,method="euclidean",diag=T,upper=T)    # distance matrix of time steps
  dmat<- as.matrix(dmat)
  size<- rep(0,ntim)
  for (i in 2:ntim) size[i]<- dmat[i-1,i]                    # taking elements below the diagonal of dmat
  size[1]<- size[2] 
  plot(time[2:ntim],size[2:ntim],xlab="Year BP",ylab="Rate of change",type="n",cex.axis=1.0,cex.lab=1.0,xaxt="s",tcl=-0.3,mgp=c(2.0,0.5,0))
  lines(time[2:ntim],size[2:ntim],lwd=0.5)  
  legend("top",expression(paste(italic("x' = "),italic(x)^0.05)),bty="n",cex=1.2)
legend("topright","(c)",cex=1.6,bty="n")
}
