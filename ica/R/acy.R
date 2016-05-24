acy <- 
  function(X,Y){
    ###### Amari-Cichocki-Yang Error
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: August 23, 2015
    
    X <- as.matrix(X)
    Y <- as.matrix(Y)
    xd <- dim(X)
    if(xd[1] != xd[2]){stop("X must be square matrix.")}
    yd <- dim(Y)
    if(yd[1] != yd[2]){stop("Y must be square matrix.")}
    if(xd[1] != yd[1]){stop("X and Y must be same dimension.")}
    A <- abs(solve(Y)%*%X)
    rowprt <- sum( (rowSums(A)/apply(A,1,max)) - 1 )
    colprt <- sum( (colSums(A)/apply(A,2,max)) - 1 )
    (rowprt + colprt) / (2*xd[1])
    
  }