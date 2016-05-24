fitted.parafac2 <-
  function(object,...){
    # Calculates Fitted Values (lists) for fit Parafac2 Models
    # Nathaniel E. Helwig (helwig@umn.edu)
    # last updated: April 9, 2015
    
    mydim <- c(NA,nrow(object$B),nrow(object$C))
    nf <- ncol(object$B)
    if(is.null(object$D)){
      fit <- vector("list",mydim[3])
      for(k in 1:mydim[3]){
        fit[[k]] <- tcrossprod(object$A$H[[k]]%*%object$A$G%*%(diag(nf)*object$C[k,]),object$B)
      }
    } else {
      nk <- sapply(object$A[[1]],nrow)
      mydim <- c(mydim,nrow(object$D))
      fit <- vector("list",mydim[4])
      for(k in 1:mydim[4]){
        fit[[k]] <- array(tcrossprod(object$A$H[[k]]%*%object$A$G%*%(diag(nf)*object$D[k,]),
                                     krprod(object$C,object$B)),dim=c(nk[k],mydim[2],mydim[3]))
      }
    }
    
    fit
    
  }