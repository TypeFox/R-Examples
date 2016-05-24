"mountains" <-
function(fcasts1,fcasts2,varnames,pts,...)
{
  # here is the setup for the 3d hills
  bw.1 <- dpill(fcasts1[,1],fcasts1[,2])
  bw.2 <- dpill(fcasts2[,1],fcasts2[,2])
             
  hill1 <- bkde2D(fcasts1, bandwidth=bw.1)
  hill2 <- bkde2D(fcasts2, bandwidth=bw.2)

  rangex <- range(hill1$x1,hill2$x1)
  rangey <- range(hill1$x2,hill2$x2)

  # here is the setup for the 2d plots at the top

  slice1.1 <- bkde(fcasts1[,1], bandwidth=bw.1)
  slice1.2 <- bkde(fcasts1[,2], bandwidth=bw.1)
  slice2.1 <- bkde(fcasts2[,1], bandwidth=bw.2)
  slice2.2 <- bkde(fcasts2[,2], bandwidth=bw.2)
  
  # plot part

   par(mai=c(0.75,0.75,0.2,0.2), cex=0.5)
   nf <- layout(matrix(c(1,3,2,4),2,2, byrow=T))
   plot(slice1.1, type="l", col=c("black"),
       xlim=range(slice1.1$x,slice2.1$x,pts[1]),
       ylim=range(slice1.1$y,slice2.1$y),
       xlab=varnames[1], ylab="Density", lwd=2)
  lines(slice2.1, col=c("red"), type="l", lwd=2)
  abline(v=pts[1])

  plot(slice1.2, type="l", col=c("black"),
       xlim=range(slice1.2$x,slice2.2$x,pts[2]),
       ylim=range(slice1.2$y,slice2.2$y),
       xlab=varnames[2], ylab="Density", lwd=2)
  lines(slice2.2, col=c("red"), type="l", lwd=2)
  abline(v=pts[2])
  
  contour(hill1$x1,hill1$x2,hill1$fhat, col=c("black"),
       xlim=range(rangex,pts[1]), ylim=range(rangey,pts[2]),
          drawlabels=F)
  
  par(new=T)
  contour(hill2$x1,hill2$x2,hill2$fhat, col=c("red"), xlim=range(rangex,pts[1]),
          ylim=range(rangey,pts[2]), xlab=varnames[1],
          ylab=varnames[2],
          drawlabels=F)
  
  points(pts[1],pts[2], pch=19)
  points(pts[3],pts[4], pch=8)
  
  par(mai=c(0.5,0.25,0.1,0.25), cex=0.4)
  persp(hill1$x1,hill1$x2,hill1$fhat,
      xlim=rangex, ylim=rangey, zlim=range(hill1$fhat,hill2$fhat),
        shade=0.3,
      xlab="", ylab="", zlab="",...)
  par(new=TRUE)
  persp(hill2$x1,hill2$x2,hill2$fhat,
        xlim=rangex, ylim=rangey, zlim=range(hill1$fhat,hill2$fhat),
        ticktype="detailed", border=c("red"), shade=0.3, cex.axis=0.5,
        xlab=varnames[1], ylab=varnames[2], zlab="Density",...)
}

