summary.bigssg <- 
  function(object,fitresid=TRUE,chunksize=10000,...){
    
    ndpts <- as.integer(object$ndf[1])
    if(fitresid){
      if(is.na(object$modelspec$rparm[1])){
        yhat <- list(fitted.values=object$fitted.values,
                     linear.predictors=object$linear.predictors)
      } else {
        chunksize <- as.integer(chunksize[1])
        if(chunksize<1L){stop("Input 'chunksize' must be positive integer.")}
        if(chunksize>=ndpts){
          yhat <- predict.bigssg(object)
        } else {
          xseq <- seq.int(1L,ndpts,by=chunksize)
          lenx <- length(xseq)
          if(xseq[lenx]<ndpts) {
            xseq <- c(xseq,ndpts+1)
            lenx <- lenx+1
          } else {
            xseq[lenx] <- xseq[lenx]+1
          }
          yhat <- vector("list",2)
          names(yhat) <- c("fitted.values","linear.predictors")
          nxvar <- length(object$xvars)
          newdata <- vector("list",nxvar)
          names(newdata) <- names(object$xvars)
          for(mm in 1:(lenx-1)){
            chunkidx <- (xseq[mm]:(xseq[mm+1]-1))
            for(k in 1:nxvar){newdata[[k]] <- object$xvars[[k]][chunkidx,]}
            ynew <- predict.bigssg(object,newdata)
            yhat$fitted.values <- c(yhat$fitted.values,ynew$fitted.values)
            yhat$linear.predictors <- c(yhat$linear.predictors,ynew$linear.predictors)
          } # end for(mm in 1:(lenx-1))
        } # end if(chunksize>=ndpts)
      } # end if(is.na(object$modelspec$rparm[1]))
      if(object$family=="binomial"){
        dev1prt <- object$yvar*log(object$yvar/(yhat$fitted.values*object$modelspec$weights))
        dev1prt[is.na(dev1prt)] <- 0
        dev2prt <- (object$modelspec$weights-object$yvar)*log((object$modelspec$weights-object$yvar)/(object$modelspec$weights-yhat$fitted.values*object$modelspec$weights))
        dev2prt[is.na(dev2prt)] <- 0
        devresid <- 2*(dev1prt+dev2prt)
        deviance <- sum(devresid,na.rm=TRUE)
        devresid <- sign(object$yvar-yhat$fitted.values*object$modelspec$weights)*sqrt(devresid)
        devresid[is.na(devresid)] <- 0
      } else if(object$family=="poisson"){
        dev1prt <- object$yvar*log(object$yvar/yhat$fitted.values)
        dev1prt[is.na(dev1prt)] <- 0
        dev2prt <- (object$yvar-yhat$fitted.values)
        devresid <- 2*(dev1prt-dev2prt)
        deviance <- sum(devresid,na.rm=TRUE)
        devresid <- sign(object$yvar-yhat$fitted.values)*sqrt(devresid)
        devresid[is.na(devresid)] <- 0
      } else if(object$family=="Gamma"){
        devresid <- 2*(-log(object$yvar/yhat$fitted.values) + (object$yvar-yhat$fitted.values)/yhat$fitted.values)
        deviance <- sum(devresid,na.rm=TRUE)
        devresid <- sign(object$yvar-yhat$fitted.values)*sqrt(devresid)
        devresid[is.na(devresid)] <- 0
      } else if(object$family=="inverse.gaussian"){
        devresid <- ((object$yvar-yhat$fitted.values)^2)/((yhat$fitted.values^2)*object$yvar)
        deviance <- sum(devresid,na.rm=TRUE)
        devresid <- sign(object$yvar-yhat$fitted.values)*sqrt(devresid)
        devresid[is.na(devresid)] <- 0
      } else if(object$family=="negbin"){
        dev1prt <- object$yvar*log(object$yvar/yhat$fitted)
        dev1prt[is.na(dev1prt)] <- 0
        dev2prt <- (object$yvar+1/object$dispersion)*log((object$yvar+1/object$dispersion)/(yhat$fitted+1/object$dispersion))
        devresid <- 2*(dev1prt-dev2prt)
        deviance <- sum(devresid,na.rm=TRUE)
        devresid <- sign(object$yvar-yhat$fitted.values)*sqrt(devresid)
        devresid[is.na(devresid)] <- 0
      }
    } else{
      yhat <- list(NULL)
      devresid <- NULL
      deviance <- NULL
    }
    sumssg <- list(call=object$call,type=object$type,fitted.values=yhat$fitted.values,
                   linear.predictors=yhat$linear.predictors,residuals=devresid,
                   deviance=deviance,dispersion=object$dispersion,n=object$ndf[1],
                   df=object$ndf[2],info=object$info,converged=object$converged,iter=object$modelspec$iter,
                   rparm=object$modelspec$rparm,lambda=object$modelspec$lambda,gammas=object$modelspec$gammas,
                   family=object$family,gcvtype=object$modelspec$gcvtype)
    class(sumssg) <- "summary.bigssg"
    return(sumssg)
    
  }