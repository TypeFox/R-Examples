larg.ellip <-
function(LSL , USL, n=25, box=FALSE, add=TRUE, xlim=xlim, ylim=ylim,zlim=zlim,

 xlab="xlab", ylab="ylab", zlab="zlab",col=2,alpha = 0.2,...){
 rx=(USL[1]-LSL[1])/2
 ry=(USL[2]-LSL[2])/2
 rz=(USL[3]-LSL[3])/2
 
    degvec <- seq(0,2*pi,length=n)
   ecoord2 <- function(p) {
     c(rx*cos(p[1])*sin(p[2]),ry*sin(p[1])*sin(p[2]),rz*cos(p[2])) }
   v <- apply(expand.grid(degvec,degvec),1,ecoord2)
v <- rbind(v,rep(1,ncol(v))) ## homogeneous
  e <- expand.grid(1:(n-1),1:n)
   i1 <- apply(e,1,function(z)z[1]+n*(z[2]-1))
   i2 <- i1+1
   i3 <- (i1+n-1) %% n^2 + 1
   i4 <- (i2+n-1) %% n^2 + 1
   i <- rbind(i1,i2,i4,i3)
 plot3d(v[1,i]+(USL[1]+LSL[1])/2,v[2,i]+(USL[2]+LSL[2])/2,v[3,i]+(USL[3]+LSL[3])/2,
 type="l",box=FALSE,col=col,add = add, alpha = alpha, xlab=xlab, ylab=ylab, zlab=zlab, xlim=xlim, ylim=ylim,zlim=zlim )
 }
