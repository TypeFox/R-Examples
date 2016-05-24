predict.bigspline <-
  function(object,newdata=NULL,se.fit=FALSE,
           effect=c("all","0","lin","non"),
           design=FALSE,smoothMatrix=FALSE,...) {
    ###### Predicts for class "bigspline" objects
    ###### Nathaniel E. Helwig (helwig@umn.edu)
    ###### Last modified: January 20, 2016
    
    ### check newdata
    effect <- effect[1]
    if(!any(effect==c("all","0","lin","non"))){stop("Must set 'effect' to one of four specified options.")}
    if(effect=="0"){
      yhat <- object$coef[1]
      if(se.fit){
        pse <- sqrt(sum(object$coef.csqrt[1,]^2))
        predcss <- list(fit=as.numeric(yhat),se.fit=pse)
      } else {
        predcss <- as.numeric(yhat)
      }
      return(predcss)
    }
    if(is.null(newdata)){
      newdata <- object$x
    } else {
      newdata <- as.matrix(newdata+0.0)
      if(ncol(newdata)>1){stop("Too many predictors in 'newdata' input.")}
    }
    nunewr <- nrow(newdata)
    nknots <- nrow(object$myknots)
    
    ### transform data and make RK matirces
    if(object$type=="cub"){
      newdata <- (newdata-object$xrng[1])/(object$xrng[2]-object$xrng[1])
      if(effect=="all"){
        Kmat <- cbind(1,newdata-0.5)
        Jmat <- (.Fortran("cubker",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 1:(nknots+2)
      } else if(effect=="non"){
        Kmat <- NULL
        Jmat <- (.Fortran("cubker",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 3:(nknots+2)
      } else {
        Kmat <- newdata-0.5
        Jmat <- NULL
        cidx <- 2
      }
    } else if(object$type=="per"){
      newdata <- (newdata-object$xrng[1])/(object$xrng[2]-object$xrng[1])
      if(effect=="all"){
        Kmat <- matrix(1,nunewr)
        Jmat <- (.Fortran("perker",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 1:(nknots+1)
      } else if(effect=="non"){
        Kmat <- NULL
        Jmat <- (.Fortran("perker",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 2:(nknots+1)
      } else {stop("There is no linear effect for periodic splines.")}
    } else if(object$type=="lin"){
      newdata <- (newdata-object$xrng[1])/(object$xrng[2]-object$xrng[1])
      if(effect=="all"){
        Kmat <- matrix(1,nunewr)
        Jmat <- (.Fortran("linker",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 1:(nknots+1)
      } else if(effect=="non"){
        Kmat <- NULL
        Jmat <- (.Fortran("linker",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 2:(nknots+1)
      } else {stop("There is no linear effect for linear splines.")}
    } else {
      if(effect=="all"){
        Kmat <- cbind(1,newdata)
        Jmat <- (.Fortran("cubkerz",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 1:(nknots+2)
      } else if(effect=="non"){
        Kmat <- NULL
        Jmat <- (.Fortran("cubkerz",newdata,object$myknots,nunewr,nknots,
                          matrix(0,nunewr,nknots),PACKAGE="bigsplines"))[[5]]
        cidx <- 3:(nknots+2)
      } else {
        Kmat <- newdata
        Jmat <- NULL
        cidx <- 2
      }
    } # end if(object$type=="cub")
    
    ### get yhat and standard errors
    yhat <- cbind(Kmat,Jmat)%*%object$coef[cidx]
    if(se.fit){pse <- sqrt(postvar(Kmat,Jmat,object$coef.csqrt[cidx,]))}
    
    ### collect new yhat
    if(se.fit | design | smoothMatrix){
      predcss <- list(fit=as.numeric(yhat))
      if(se.fit){predcss <- c(predcss,list(se.fit=pse))}
      if(design){predcss <- c(predcss,list(X=cbind(Kmat,Jmat),ix=cidx))}
      if(smoothMatrix){predcss <- c(predcss,list(S=tcrossprod(cbind(Kmat,Jmat)%*%object$coef.csqrt[cidx,])))}
    } else {
      predcss <- as.numeric(yhat)
    }
    return(predcss)
    
  }