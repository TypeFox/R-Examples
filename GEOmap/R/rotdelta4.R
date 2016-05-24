rotdelta4<-function(delta)
{
  rotmat=matrix(ncol=4, nrow=4)

  cdel = cos(delta*pi/180)
  sdel = sin(delta*pi/180)

  rotmat[1,1] <- cdel
  rotmat[1,2] <- sdel
  rotmat[1,3] <- 0.0
  rotmat[1,4] <- 0.0

  rotmat[2,1] <- -sdel
  rotmat[2,2] <- cdel
  rotmat[2,3] <- 0
  rotmat[2,4] <- 0.0

  rotmat[3,1] <- 0.0
  rotmat[3,2] <- 0
  rotmat[3,3] <-  1
  rotmat[3,4] <- 0.0

  rotmat[4,1] <- 0.0
  rotmat[4,2] <- 0.0
  rotmat[4,3] <- 0.0
  rotmat[4,4] <- 1.0
  
  return(rotmat)
}

