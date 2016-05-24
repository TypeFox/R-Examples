"print.maboost" <-
  function(x,...){
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
    if(!is.null(cl <- x$call)) {
      names(cl)[2] <- ""
      cat("Call:\n")
      dput(cl)
    }
    
    g=x$fit
    tab=x$confusion
    errm=1-sum(diag(tab))/length(x$actual)
    
    cat("\nLoss:", x$model$lossObj$loss,"Method:", x$model$lossObj$type,"\  Iteration:", x$iter,"\n")
    cat("\nFinal Confusion Matrix for Data:\n")
    print(tab)
    cat("\nTrain Error:", round(errm, digits=3),"\n")
    if(x$iter<=5){
      j<-which.min(x$model$oob.str$oobm.err)
    }else{
      j<-which.min(x$model$oob.str$oobm.err[-c(1:5)])+5
    }
    cat("\nOut-Of-Bag Error: ",round(x$model$oob.str$oobm.err[j],3)," iteration=",j,"\n")
    cat("\nAdditional Estimates of number of iterations:\n\n")
    errs<-as.data.frame(x$model$errs)
    k<-dim(errs)[2]/2
    names(errs)<-paste(names(errs),sort(c(1:k,1:k)),sep="")
    est.m<-apply(errs,2,which.min)
    print(est.m)
    cat("\n")
    invisible(x)
  }
