predict.bigtps <-
  function(object,newdata=NULL,se.fit=FALSE,
           effect=c("all","0","lin","non"),
           design=FALSE,smoothMatrix=FALSE,...) {
    ###### Predicts for class "bigtps" objects
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: March 10, 2015    
    
    ### check newdata
    effect <- effect[1]
    nx <- ncol(object$x)
    if(any(effect==c("all","0","lin","non"))==FALSE){
      stop("Must set 'effect' to one of four specified options.")
    }
    if(effect=="0"){
      yhat <- object$coef[1]
      if(se.fit){
        pse <- sqrt(sum(object$coef.csqrt[1,]^2))
        predtps <- list(fit=as.numeric(yhat),se.fit=pse)
      } else {
        predtps <- as.numeric(yhat)
      }
      return(predtps)
    }
    if(is.null(newdata)){
      newdata <- object$x
    } else {
      newdata <- as.matrix(newdata+0.0)
      if(ncol(newdata)!=nx){stop("Wrong number of predictors in 'newdata' input.")}
    }
    nunewr <- nrow(newdata)
    nknots <- nrow(object$myknots)
    
    ### make RK matirces
    if(nx<3){gconst=1} else{gconst=-1}
    if(effect=="all"){
      Kmat <- cbind(1,newdata)
      Jmat <- gconst*(.Fortran("tpsker",newdata,object$myknots,nunewr,nx,nknots,matrix(0,nunewr,nknots)))[[6]]
      cidx <- 1:(nknots+nx+1)
    } else if(effect=="lin"){
      Kmat <- newdata
      Jmat <- NULL
      cidx <- 2:(nx+1)
    } else {
      Kmat <- NULL
      Jmat <- gconst*(.Fortran("tpsker",newdata,object$myknots,nunewr,nx,nknots,matrix(0,nunewr,nknots)))[[6]]
      cidx <- (nx+2):(nknots+nx+1)
    } # end if(effect=="all")
    
    ### get yhat and standard errors
    yhat <- cbind(Kmat,Jmat)%*%object$coef[cidx]
    if(se.fit){pse <- sqrt(postvar(Kmat,Jmat,object$coef.csqrt[cidx,]))}
    
    ### collect new yhat
    if(se.fit | design | smoothMatrix){
      predtps <- list(fit=as.numeric(yhat))
      if(se.fit){predtps <- c(predtps,list(se.fit=pse))}
      if(design){predtps <- c(predtps,list(X=cbind(Kmat,Jmat),ix=cidx))}
      if(smoothMatrix){predtps <- c(predtps,list(S=tcrossprod(cbind(Kmat,Jmat)%*%object$coef.csqrt[cidx,])))}
    } else {
      predtps <- as.numeric(yhat)
    }
    return(predtps)
    
  }