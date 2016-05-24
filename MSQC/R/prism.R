prism <-
function(LSL=LSL,USL=USL, add=TRUE,xlim=xlim, ylim=ylim,zlim=zlim,...){ 

 xf<-c(LSL[1],LSL[1],LSL[1],LSL[1],LSL[1])
 yf<-c(LSL[2],USL[2],USL[2],LSL[2],LSL[2])
 zf<-c(LSL[3],LSL[3],USL[3],USL[3],LSL[3])
 
 plot3d(xf,yf,zf, type="l",add=add,...)
 
  xf<-c(USL[1],USL[1],USL[1],USL[1],USL[1])
 yf<-c(LSL[2],USL[2],USL[2],LSL[2],LSL[2])
 zf<-c(LSL[3],LSL[3],USL[3],USL[3],LSL[3])
 
 plot3d(xf,yf,zf, type="l",add=add,...)
 
 plot3d(c(LSL[1], USL[1]), LSL[2], LSL[3], type="l",add=add, ...)
 plot3d(c(LSL[1], USL[1]), USL[2], LSL[3], type="l",add=add, ...)
 plot3d(c(LSL[1], USL[1]), LSL[2], USL[3], type="l",add=add, ...)
 plot3d(c(LSL[1], USL[1]), USL[2], USL[3], type="l",add=add, ...)
 
 }
