reorder.tucker <-
  function(x, neworder, mode="A", ...){
    # Reorder Factors of fit Tucker model
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: August 19, 2015    
    
    # check mode
    mode <- mode[1]
    if(is.null(x$D)){
      if(!any(mode==c("A","B","C"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', or 'C' for 3-way Tucker")
    } else {
      if(!any(mode==c("A","B","C","D"))) stop("Incorrect input for 'mode'. Must set to 'A', 'B', 'C', or 'D' for 4-way Tucker")
    }
    
    # check neworder and reorder factors
    neworder <- as.integer(neworder)
    if(mode=="A"){
      
      nfac <- ncol(x$A)
      if(length(neworder)!=nfac) stop(paste("Incorrect input for 'neworder'. Must have length equal to number of factors for Mode",mode))
      if(!identical(seq(1L,nfac),sort(neworder))) stop(paste("Incorrect input for 'neworder'. Must be unique integers in range of 1 to",nfac))
      x$A <- x$A[,neworder]
      if(is.null(x$D)) {
        x$G <- x$G[neworder,,]
      } else {
        x$G <- x$G[neworder,,,]
      }
      
    } else if(mode=="B"){
      
      nfac <- ncol(x$B)
      if(length(neworder)!=nfac) stop(paste("Incorrect input for 'neworder'. Must have length equal to number of factors for Mode",mode))
      if(!identical(seq(1L,nfac),sort(neworder))) stop(paste("Incorrect input for 'neworder'. Must be unique integers in range of 1 to",nfac))
      x$B <- x$B[,neworder]
      if(is.null(x$D)) {
        x$G <- x$G[,neworder,]
      } else {
        x$G <- x$G[,neworder,,]
      }
      
    } else if(mode=="C"){
      
      nfac <- ncol(x$C)
      if(length(neworder)!=nfac) stop(paste("Incorrect input for 'neworder'. Must have length equal to number of factors for Mode",mode))
      if(!identical(seq(1L,nfac),sort(neworder))) stop(paste("Incorrect input for 'neworder'. Must be unique integers in range of 1 to",nfac))
      x$C <- x$C[,neworder]
      if(is.null(x$D)) {
        x$G <- x$G[,,neworder]
      } else {
        x$G <- x$G[,,neworder,]
      }
      
    } else if(mode=="D"){
      
      nfac <- ncol(x$D)
      if(length(neworder)!=nfac) stop(paste("Incorrect input for 'neworder'. Must have length equal to number of factors for Mode",mode))
      if(!identical(seq(1L,nfac),sort(neworder))) stop(paste("Incorrect input for 'neworder'. Must be unique integers in range of 1 to",nfac))
      x$D <- x$D[,neworder]
      x$G <- x$G[,,,neworder]
      
    }
    
    return(x)
    
  }