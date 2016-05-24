fitted.parafac <-
  function(object,...){
    # Calculates Fitted Values (arrays) for fit Parafac Models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    mydim <- c(nrow(object$A),nrow(object$B),nrow(object$C))
    if(is.null(object$D)){
      fit <- array(tcrossprod(object$A,krprod(object$C,object$B)),dim=mydim)
    } else {
      mydim <- c(mydim,nrow(object$D))
      fit <- array(tcrossprod(object$A,krprod(object$D,krprod(object$C,object$B))),dim=mydim)
    }
    
    fit
    
  }