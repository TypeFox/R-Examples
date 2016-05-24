resign.indscal <-
  function(x, mode="B", newsign=1, ...){
    # Resigns Weights of fit INDSCAL model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: August 19, 2015    
    
    # check newsign
    nfac <- ncol(x$B)
    newsign <- sign(newsign)
    if(length(newsign)!=nfac) newsign <- rep(newsign[1],nfac)
    
    # resign factors
    Bsign <- sign(colMeans(x$B^3))
    svec <- newsign*Bsign
    if(nfac==1L) { Smat <- matrix(svec) } else { Smat <- diag(svec) }
    x$B <- x$B %*% Smat
    return(x)
    
  }