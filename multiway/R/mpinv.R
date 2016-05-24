mpinv <- 
  function(X,tol=.Machine$double.eps){
    # Moore-Penrose Pseudoinverse
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: October 10, 2015
    
    X <- as.matrix(X)
    xsvd <- svd(X)
    nze <- sum( xsvd$d > (tol*xsvd$d[1]) )
    if(nze > 1L){
      return( xsvd$v[,1:nze] %*% diag(1/xsvd$d[1:nze]) %*% t(xsvd$u[,1:nze]) )
    } else {
      return( outer(xsvd$v[,1],xsvd$u[,1]) / xsvd$d[1] )
    }
    
}