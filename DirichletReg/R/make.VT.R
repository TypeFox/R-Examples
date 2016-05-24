D2R <- function(degrees) pi*degrees/180



t3d <- function(XYZ, ViewTrans3D) trans3d(XYZ[,1], XYZ[,2], XYZ[,3], ViewTrans3D)



make.VT <- function(theta=0,
                    phi=15,
                    r=sqrt(3),
                    d=1,
                    expand=1,
                    origin=c(0,0,0),
                    scale=c(1,1,1)){

  xc <- origin[1]
  yc <- origin[2]
  zc <- origin[3]

  xs <- scale[1]
  ys <- scale[2]
  zs <- scale[3]

  VT <- diag(4)   # initialize

  TT <- diag(4)   # center @ origin
  TT[4,1:3] <- c(-xc, -yc, -zc)
  VT <- VT %*% TT


  TT <- diag(4)   # scale extents to [-1, 1]
  TT[cbind(1:3,1:3)] <- c(1/xs, 1/ys, expand/zs)
  VT <- VT %*% TT


  TT <- diag(4)   # rotate x-y plane to horizontal
  TT[cbind(c(2,3,3,2),c(2,2,3,3))] <- c(cos(D2R(-90)), -sin(D2R(-90)),
                                        cos(D2R(-90)),  sin(D2R(-90)))
  VT <- VT %*% TT


  TT <- diag(4)   # azimuthal rotation (theta)
  TT[cbind(c(1,3,3,1),c(1,1,3,3))] <- c(cos(D2R(-theta)),  sin(D2R(-theta)),
                                        cos(D2R(-theta)), -sin(D2R(-theta)))
  VT <- VT %*% TT


  TT <- diag(4)   # elevation rotation (phi)
  TT[cbind(c(2,3,3,2),c(2,2,3,3))] <- c(cos(D2R(phi)), -sin(D2R(phi)),
                                        cos(D2R(phi)),  sin(D2R(phi)))
  VT <- VT %*% TT


  TT <- diag(4)   # translate eyepoint to origin
  TT[4,1:3] <- c(0, 0, -r-d)
  VT <- VT %*% TT


  TT <- diag(4)   # perspective
  TT[3,4] <- -1/d
  VT <- VT %*% TT

  return(VT)
}
