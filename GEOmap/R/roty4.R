roty4<-function(vec)
{ 
  rotmat=matrix(ncol=4, nrow=4)

  r = sqrt(vec[1]^2 +vec[2]^2 + vec[3]^2)
  
  cee = vec/r
  d = sqrt(cee[2]^2 + cee[3]^2)
  
  rotmat[1,1] <- d
  rotmat[1,2] <- 0
  rotmat[1,3] <- cee[1]
  rotmat[1,4] <- 0.0

  rotmat[2,1] <- 0.0
  rotmat[2,2] <- 1
  rotmat[2,3] <- 0
  rotmat[2,4] <- 0.0

  rotmat[3,1] <- -cee[1]
  rotmat[3,2] <- 0
  rotmat[3,3] <-  d
  rotmat[3,4] <- 0.0

  rotmat[4,1] <- 0.0
  rotmat[4,2] <- 0.0
  rotmat[4,3] <- 0.0
  rotmat[4,4] <- 1.0

  
  
  return(rotmat)
}
