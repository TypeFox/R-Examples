fitted.tucker <-
  function(object,...){
    # Calculates Fitted Values (arrays) for fit Tucker Models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    mydim <- c(nrow(object$A),nrow(object$B),nrow(object$C))
    nf <- c(ncol(object$A),ncol(object$B),ncol(object$C))
    if(is.null(object$D)){
      fit <- array(tcrossprod(object$A%*%matrix(object$G,nf[1],nf[2]*nf[3]),kronecker(object$C,object$B)),dim=mydim)
    } else {
      mydim <- c(mydim,nrow(object$D))
      nf <- c(nf,ncol(object$D))
      fit <- array(tcrossprod(object$A%*%matrix(object$G,nf[1],prod(nf[2:4])),
                              kronecker(object$D,kronecker(object$C,object$B))),dim=mydim)
    }
    
    fit
    
  }