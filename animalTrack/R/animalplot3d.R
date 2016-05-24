animalplot3d <-
function(roll, pitch, yaw, angle = "degree", add=FALSE, ...){
  
  suml <- length(roll) + length(pitch) + length(yaw)
  
  if (suml != 3) stop("roll, pitch, and yaw must be 1 element objects")
    
  # Set up the viewing frame and viewer 
  goodview <- matrix(c(0.8451259,-0.1130982,0.5224662,0,0.5344371,0.2003249,-0.8211254,0,
                       -0.0117952,0.9731796,0.2297437,0,0,0,0,1),ncol=4)
  
  # Now lets rotate the coordinate system to NED
  Yaw90ccw <- matrix(c(0,-1,0,1,0,0,0,0,1),ncol=3)  # counter clockwise yaw (rotation about the z axis)
  ReflectY <- matrix(c(1,0,0,0,-1,0,0,0,1),ncol=3)  #reflection of the y axis about the xz plane
  ReflectZ <- matrix(c(1,0,0,0,1,0,0,0,-1),ncol=3)  #reflection of the z axis about the xy plane
  staticM <- matrix(c(1,0,0,0,1,0,0,0,1),ncol=3)  #identity matrix
  staticM2 <- staticM %*% ReflectY %*% ReflectZ %*% Yaw90ccw  #rotate the staticM matrix to the animal frame of reference (i.e. 
                                                              # north-east-down NED)
  
  #degree to radians, if need be.
  if (angle == "degree"){ roll <- roll * (pi/180) ; pitch <- -(pitch) * (pi/180) ; yaw <- -(yaw) * (pi/180) }
  if (angle == "radian"){ roll <- roll ; pitch <- -(pitch) ; yaw <- -(yaw) }
  
  # Roll, pitch, and yaw rotation matrices
  R = matrix(c(1,0,0,0,cos(roll),sin(roll),0,-sin(roll),cos(roll)),ncol=3)
  P = matrix(c(cos(pitch),0,-sin(pitch),0,1,0,sin(pitch),0,cos(pitch)),ncol=3)
  Y = matrix(c(cos(yaw),sin(yaw),0,-sin(yaw),cos(yaw),0,0,0,1),ncol=3)
  Cprime <- R %*% P %*% Y %*% staticM
  Cprime2 <- Cprime %*% ReflectY %*% ReflectZ %*% Yaw90ccw
  
  # rgl plot  
  textoffset <- 1.2 ; limits <- c(-1.4,1.4) ; smallradius <- 0.02 
  plot3d(staticM2[1,],staticM2[2,],staticM2[3,],xlim=limits,ylim=limits,zlim=limits,type="s",col="red",
         xlab="Y Axis",ylab="X Axis",zlab="Z Axis", box=FALSE, radius=smallradius, add=add, ...)  
  texts3d(staticM2[1,]*textoffset,texts=c("X")) ; texts3d(staticM2[2,]*textoffset,texts=c("Y")) ;texts3d(staticM2[3,]*textoffset,texts=c("Z"))
  segments3d(cbind(c(0,staticM2[1,1]),c(0,staticM2[2,1]),c(0,staticM2[3,1])),col="red",lwd=1)
  segments3d(cbind(c(0,staticM2[1,2]),c(0,staticM2[2,2]),c(0,staticM2[3,2])),col="red",lwd=1)
  segments3d(cbind(c(0,staticM2[1,3]),c(0,staticM2[2,3]),c(0,staticM2[3,3])),col="red",lwd=1)
  segments3d(cbind(c(0,Cprime2[1,1]),c(0,Cprime2[1,2]),c(0,Cprime2[1,3])),col="blue",lwd=2) 
  segments3d(cbind(c(0,Cprime2[2,1]),c(0,Cprime2[2,2]),c(0,Cprime2[2,3])),col="blue",lwd=2)
  segments3d(cbind(c(0,Cprime2[3,1]),c(0,Cprime2[3,2]),c(0,Cprime2[3,3])),col="blue",lwd=2)
  texts3d(Cprime2[1,]*textoffset,texts=c("X"),cex=1.5) ; texts3d(Cprime2[2,]*textoffset,texts=c("Y"),cex=1.5) ;texts3d(Cprime2[3,]*textoffset,texts=c("Z"),cex=1.5)
  rgl.spheres( Cprime2[1,],col="blue",radius=0.07) ; rgl.spheres( Cprime2[2,],col="blue",radius=smallradius) ; rgl.spheres( Cprime2[3,],col="blue",radius=smallradius)
  rgl.viewpoint(userMatrix = goodview)

  # Return the rotation matrix and the homogenous rotatation matrix used to make the rgl plot
  return(list("rotmat" = Cprime2,"hrotmat" = asHomogeneous(Cprime2)))
  
}

