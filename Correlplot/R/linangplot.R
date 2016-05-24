linangplot <- function(x,y,tmx=NULL,tmy=NULL,...) {
  r <- cor(x,y)
#  cat("correlaton:",r,"\n")
  ang <- pi/2*(1-r)
#  cat("angle: ",ang,"radians \n")
  deg <- rad2degree(ang)
#  cat("angle: ",deg,"degrees \n")
  cosang <- cos(ang)
  rl <- (-2/pi)*ang+1
  b21 <- cosang
  b22 <- sqrt(1-b21^2)
  b2 <- c(b21,b22)
#  print(b2)
  b1 <- c(1,0)
  B <- cbind(b1,b2)
  Binv <- ginv(B)

  xs <- scale(x)
  ys <- scale(y)
  Ft <- cbind(xs,ys)%*%Binv
  opar <- par("yaxt"="n","xaxt"="n")
  plot(Ft[,1],Ft[,2],pch=19,xlab="",ylab="",asp=1,...)
#  textxy(Ft[,1],Ft[,2],1:nrow(Ft))
#  arrows(0,0,B[1,1],B[2,1],length=0.1,lwd=2)
#  arrows(0,0,B[1,2],B[2,2],length=0.1,lwd=2)

  if(is.null(tmx)) {
    tms <- seq(-5,5,1)
    Calibrate.x <- calibrate(B[,1],xs,tms,tmlab=tms,Ft,axislab="",graphics=TRUE,axiscol="black",where=3,verb=FALSE)
  }
  else {
    tmc <- tmx - mean(x)
    tms <- tmc/sqrt(var(x))
    Calibrate.x <- calibrate(B[,1],xs,tms,tmlab=tmx,Ft,axislab="",graphics=TRUE,axiscol="black",where=3,verb=FALSE)
  }
  if(is.null(tmy)) {
    tms <- seq(-5,5,1)
    Calibrate.y <- calibrate(B[,2],ys,tms,tmlab=tms,Ft,axislab="",graphics=TRUE,axiscol="black",where=2,verb=FALSE)
  }
  else {
     tmc <- tmy - mean(y)
     tms <- tmc/sqrt(var(y))
     Calibrate.y <- calibrate(B[,2],ys,tms,tmlab=tmy,Ft,axislab="",graphics=TRUE,axiscol="black",where=2,verb=FALSE)
  }
  par(opar)
  return(list(Xt=Ft,B=B,r=r,angledegrees=deg,angleradians=ang))
}

