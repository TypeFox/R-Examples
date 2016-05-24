varmx <- function(amat, normalize=FALSE) {

  #  Does a VARIMAX rotation of a principal components solution

  #  Arguments:
  #  AMAT      ...  N by K matrix of harmonic values
  #  NORMALIZE ... either TRUE or FALSE.  If TRUE, the columns of AMAT
  #                are normalized prior to computing the rotation 
  #                matrix.  However, this is seldom needed for 
  #                functional data.

  #  Returns:
  #  ROTM  ...  Rotation matrixed loadings

  #  Last modified 22 October by Jim Ramsay

  n    <- nrow(amat)
  k    <- ncol(amat)
  rotm <- diag(k)
  onek <- matrix(1,1,k)

  if (k == 1) return(rotm)

  #  normalize loadings matrix

  if (normalize) {
      hvec <- as.matrix(apply(amat^2, 1, var))
      amat <- amat/(sqrt(hvec) %*% onek)
  }

  eps  <- 0.0011
  ccns <- 0.7071068

  varold <- 0
  varnow <- sum(apply(amat^2, 2, var))

  iter <- 0
  while (abs(varnow - varold) > 1e-7 && iter <= 50) {
    iter  <- iter + 1
    for (j in 1:(k-1)) for (l in (j+1):k) {
      avecj  <- amat[,j]
      avecl  <- amat[,l]
      uvec   <- avecj^2 - avecl^2
      tvec   <- 2*avecj*avecl
      aa <- sum(uvec)
      bb <- sum(tvec)
      cc <- sum(uvec^2 - tvec^2)
      dd <- 2*sum(uvec*tvec)
      tval <- dd - 2*aa*bb/n
      bval <- cc - (aa^2 - bb^2)/n

      if (tval == bval) {
        sin4t <- ccns
        cos4t <- ccns
      }

      if (tval < bval) {
        tan4t <- abs(tval/bval)
        if (tan4t >= eps) {
          cos4t <- 1/sqrt(1+tan4t^2)
          sin4t <- tan4t*cos4t
        } else {
          if (bval < 0) {
            sin4t <- ccns
            cos4t <- ccns
          } else {
            sin4t <- 0
            cos4t <- 1
          }
        }
      }

      if (tval > bval) {
        ctn4t <- abs(tval/bval)
        if (ctn4t >= eps) {
          sin4t <- 1/sqrt(1+ctn4t^2)
          cos4t <- ctn4t*sin4t
        } else {
          sin4t <- 1
          cos4t <- 0
        }
      }

      cos2t <- sqrt((1+cos4t)/2)
      sin2t <- sin4t/(2*cos2t)
      cost  <- sqrt((1+cos2t)/2)
      sint  <- sin2t/(2*cost)
      if (bval > 0) {
        cosp <- cost
        sinp <- sint
      } else {
        cosp <- ccns*(cost + sint)
        sinp <- ccns*abs(cost - sint)
      }
      if (tval <= 0) sinp <- -sinp

      amat[,j] <-  avecj*cosp + avecl*sinp
      amat[,l] <- -avecj*sinp + avecl*cosp
      rvecj    <- rotm[,j]
      rvecl    <- rotm[,l]
      rotm[,j] <-  rvecj * cosp + rvecl * sinp
      rotm[,l] <- -rvecj * sinp + rvecl * cosp

    }
    varold <- varnow
    varnow <- sum(apply(amat^2,2,var))
  }

  return( rotm )
}
