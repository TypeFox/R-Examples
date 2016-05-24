"mean.prodlim" <- function(x,
                           times,
                           newdata,
                           ...){

  if (!(x$model %in% c("survival","competing.risks"))) stop("no mean(.prodlim) method available for this object.")
  if(x$covariate.type==1) stop("No covariates for computing mean survival.") 
  
  jump.times <- x$time
  if (missing(times)) times <- x$time
  times <- sort(unique(times))
  ntimes <- length(times)
  if (missing(newdata)) newdata <- eval(x$call$data)
  surv.frame <- predict(x,newdata=newdata,time=times,level.chaos=1,mode="matrix",type="surv")
  smean <- apply(surv.frame,2,mean,na.rm=TRUE)
  marginal.fit <- prodlim(update.formula(formula(x$formula),"~1"),data=x$data)
  out <- marginal.fit
  out$surv <- smean
  out$covariate.type <- 1
  class(out) <- c("prodlim","mean")
  out
}

