ed2sp <-
  function(D,checks=TRUE){
    # Euclidean Distance to Scalar Product
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    # initial info
    D <- as.matrix(D)
    n <- nrow(D)
    if(checks){
      if(n!=ncol(D)){stop("D must be square matrix.")}
      if(mean(abs(diag(D)))>.Machine$double.eps){stop("D must have zeros on diagonal.")}
      if(!isSymmetric(D)){stop("D must be symmetric matrix.")}
    }
    
    # transform ED to SP
    Dsq <- D^2
    a <- rowMeans(Dsq)
    b <- colMeans(Dsq)
    c <- mean(Dsq)
    (-0.5)*( Dsq - matrix(1,n,1)%*%a - b%*%matrix(1,1,n) + c )
    
  }