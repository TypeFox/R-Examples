trans4<-function(vec)
{
                                        #   translation transormation
  rotmat=matrix(ncol=4, nrow=4)

  rotmat[1,1] <- 1
  rotmat[1,2] <- 0
  rotmat[1,3] <- 0
  rotmat[1,4] <- 0

  rotmat[2,1] <- 0
  rotmat[2,2] <- 1
  rotmat[2,3] <- 0
  rotmat[2,4] <- 0

  rotmat[3,1] <- 0
  rotmat[3,2] <- 0
  rotmat[3,3] <- 1
  rotmat[3,4] <- 0

  rotmat[4,1] <- vec[1]
  rotmat[4,2] <- vec[2]
  rotmat[4,3] <- vec[3]
  rotmat[4,4] <- 1.0
  return(rotmat)
}
#############################
