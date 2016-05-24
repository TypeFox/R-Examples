"survest.aalen" <- function(fit,newdata,times,...){
    ##  Time-varying hazard 
    time.coef <- data.frame(fit$cum)
    ntime <- nrow(time.coef)
    fittime <- time.coef[,1,drop=TRUE]
    ntimevars <- ncol(time.coef)-2
    time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
    nobs <- nrow(newdata)
    hazard <- .C("survest_cox_aalen",
                 timehazard=double(ntime*nobs),
                 as.double(unlist(time.coef[,-1])),
                 as.double(unlist(time.vars)),
                 as.integer(ntimevars+1),
                 as.integer(nobs),
                 as.integer(ntime),package="pecDev")$timehazard
    hazard <- matrix(hazard,ncol=ntime,nrow=nobs,dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
    surv <- pmin(exp(-hazard),1)
    if (missing(times)) times <- sort(unique(fittime))
    pred <- surv[,prodlim::sindex(jump.times=fittime,eval.times=times)]
    #    class(pred) <- c("survest","cox.aalen")
    pred
  }
