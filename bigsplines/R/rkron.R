rkron <-
function(X,Y){
  # row-wise Kronecker product
  
  if(is.null(X) | is.null(Y)){
    Z=NULL
  } else {
    X <- as.matrix(X)
    xd <- dim(X)
    Y <- as.matrix(Y)
    yd <- dim(Y)
    if(xd[1]!=yd[1]){stop("X and Y must have the same number of rows.")}
    Z <- matrix(0,xd[1],xd[2]*yd[2])
    for(jj in 1:yd[2]){
      zind <- (1+(jj-1)*xd[2]):(xd[2]+(jj-1)*xd[2])
      Z[,zind] <- X*Y[,jj]
    }
  }
  Z
  
}
