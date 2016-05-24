"survest.cox.aalen" <- function(fit,newdata,times,...){
  
  ##  The time-constant effects first
  const <- c(fit$gamma)
  names(const) <- substr(dimnames(fit$gamma)[[1]],6,nchar(dimnames(fit$gamma)[[1]])-1)
  constant.part <- t(newdata[,names(const)])*const
  constant.part <- exp(colSums(constant.part))
  
  ##  Then extract the time-varying effects
  time.coef <- data.frame(fit$cum)
  ntime <- nrow(time.coef)
  fittime <- time.coef[,1,drop=TRUE]
  ntimevars <- ncol(time.coef)-2
  
  time.vars <- cbind(1,newdata[,names(time.coef)[-(1:2)],drop=FALSE])
  nobs <- nrow(newdata)
  
  time.part <- .C("survest_cox_aalen",timehazard=double(ntime*nobs),as.double(unlist(time.coef[,-1])),as.double(unlist(time.vars)),as.integer(ntimevars+1),as.integer(nobs),as.integer(ntime),package="pecDev")$timehazard
  
  time.part <- matrix(time.part,
                      ncol=ntime,
                      nrow=nobs,
                      dimnames=list(1:nobs,paste("TP",1:ntime,sep="")))
  
  surv <- pmin(exp(-time.part*constant.part),1)
  
  if (missing(times)) times <- sort(unique(fittime))
  pred <- surv[,prodlim::sindex(fittime,times)]
  class(pred) <- c("survest","cox.aalen")
  pred
}
