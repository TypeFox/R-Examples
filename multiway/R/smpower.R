smpower <- 
  function(X,power=0.5,tol=.Machine$double.eps){
    # Symmetric Matrix Power
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    Xeig <- eigen(X,symmetric=TRUE)
    absval <- abs(Xeig$val)
    nze <- which( absval > (tol*max(absval)) )
    Xpwr <- tcrossprod(Xeig$vec[,nze]%*%(diag(length(nze))*(Xeig$val[nze]^power)),Xeig$vec[,nze])
    Xpwr
    
  }