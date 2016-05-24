"summary.maboost" <-
  function(object,n.iter=NULL,...){
    if(!inherits(object,"maboost")){
      stop("Error:  Object is not of class maboost")
    }
    kapstat<-function(tab=diag(2) ){
      if(dim(tab)[1]==1){
        return(0)
      }
      if(dim(tab)[1]==dim(tab)[2]){
        rs<-apply(tab,2,sum)
        cs<-apply(tab,1,sum)
        N<-sum(rs)
        E<-sum(rs*cs)/N^2
        O<-sum(diag(tab))/N
        return( (O-E)/(1-E) )
      }else{
        return(NA)
      }
    }
    if(!is.null(cl <- object$call)) {
      names(cl)[2] <- ""
      cat("Call:\n")
      dput(cl)
    }
    iter=object$iter
    if(!is.null(n.iter)){
      if(n.iter<iter)
        iter=n.iter
    }
    g=object$fit
    tab=object$confusion
    errm=object$model$errs[iter,]
    cat("\nLoss:",object$model$lossObj$loss,"Method:", object$model$lossObj$type,"\  Iteration:", iter,"\n")
    
    cat("\nTraining Results\n")
    
    cat("\nAccuracy:", round(1-errm[1],digits=3),
        "Kappa:",round(1-errm[2], digits=3) ,"\n\n")
    k<-length(errm)/2
    if(k>1){
      l<-3
      cat("Testing Results\n")
      for(i in 2:k){
        cat("\nAccuracy:", round(1-errm[l], digits=3))
        l<-l+1
        cat(" Kappa:", round(1-errm[l], digits=3),"\n")
        l<-l+1
      }
      cat("\n\n")
    }
  }
