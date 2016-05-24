summary.bigtps <- 
  function(object,fitresid=TRUE,chunksize=10000,...){
    
    ndpts <- as.integer(object$ndf[1])
    if(fitresid){
      if(is.na(object$rparm[1])){
        yhat <- object$fitted.values
      } else {
        chunksize <- as.integer(chunksize[1])
        if(chunksize<1L){stop("Input 'chunksize' must be positive integer.")}
        if(chunksize>=ndpts){
          yhat <- predict.bigtps(object)
        } else {
          xseq <- seq.int(1L,ndpts,by=chunksize)
          lenx <- length(xseq)
          if(xseq[lenx]<ndpts) {
            xseq <- c(xseq,ndpts+1)
            lenx <- lenx+1
          } else {
            xseq[lenx] <- xseq[lenx]+1
          }
          yhat <- NULL
          for(mm in 1:(lenx-1)){
            chunkidx <- (xseq[mm]:(xseq[mm+1]-1))
            yhat <- c(yhat,predict.bigtps(object,newdata=object$x[chunkidx,]))
          } # end for(mm in 1:(lenx-1))
        } # end if(chunksize>=ndpts)
      } # end if(is.na(object$rparm[1]))
      resid <- object$y-yhat
    } else{
      yhat <- resid <- NULL
    }
    sumtps <- list(call=NA,type="tps",fitted.values=yhat,residuals=resid,
                   sigma=object$sigma,n=ndpts,df=object$ndf[2],info=object$info,
                   converged=NA,iter=NA,rparm=object$rparm,lambda=object$lambda)
    class(sumtps) <- "summary.bigtps"
    return(sumtps)
    
  }