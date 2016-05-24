FLLat <- function(Y,J=min(15,floor(ncol(Y)/2)),B="pc",lam1,lam2,
                  thresh=10^(-4),maxiter=100,maxiter.B=1,maxiter.T=1) {
  
  ## Error checking parameters.
  CheckPars(Y,J,B,lam1,lam2,thresh,maxiter,maxiter.B,maxiter.T)

  ## Setting weight constraint.
  sT <- 1
  
  ## Initializing Beta and Theta.
  if (is.matrix(B)) {
    old.B <- B
  } else if (B=="pc") {
    Y.cen <- scale(Y,scale=FALSE)
    Y.v <- svd(Y.cen)$v
    old.B <- Y.cen%*%(Y.v[,1:J,drop=FALSE])
  } else if (B=="rand") {
    old.B <- Y[,sample(ncol(Y),J),drop=F]
  }
  old.T <- matrix(0,nrow=J,ncol=ncol(Y))

  ## Running FLLat.
  result <- .Call(LatL2C,Y,as.integer(J),old.B,old.T,as.double(lam1),
                  as.double(lam2),as.double(thresh),as.integer(maxiter),
                  as.integer(maxiter.B),as.integer(maxiter.T),
                  as.double(sT))
  return(result)

}
