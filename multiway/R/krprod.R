krprod <- 
  function(X,Y){
    # Khatri-Rao product (columnwise Kronecker product)
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    xdim <- dim(X)
    ydim <- dim(Y)
    if(xdim[2]!=ydim[2]){stop("X and Y must have same number of columns.")}
    XkrY <- matrix(0,ydim[1]*xdim[1],xdim[2])
    for(u in 1:xdim[2]){
      XkrY[,u] <- kronecker(X[,u],Y[,u])
    }
    XkrY
  }