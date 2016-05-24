###########################################################
##### cross-validation for log penalized regression  ######
###########################################################

## just an R loop that calls gamlr
cv.gamlr <- function(x, y, nfold=5, foldid=NULL, verb=FALSE, cl=NULL, ...){
  
  full <- gamlr(x,y, ...)
  fam <- full$family
  y <- checky(y,fam)

  nobs <- full$nobs
  if(is.null(foldid)){
    nfold <- min(nfold,nobs)
    foldsize <- ceiling(nobs/nfold)
    foldid <- rep.int(1:nfold,times=foldsize)[sample.int(nobs)]
  } else  stopifnot(length(foldid)==nobs)
  foldid <- factor(foldid)
  nfold <- nlevels(foldid)

  argl <- list(...)
  if(!is.null(argl$shift)) shift <- argl$shift
  lambda <- as.double(full$lambda)
  argl$lambda.start <- lambda[1]
  argl$nlambda <- length(lambda)
  argl$lambda.min.ratio <- tail(lambda,1)/lambda[1]

  ## remove any pre-calculated summaries
  argl$vxx <- argl$vxsum <- argl$vxy <- argl$xbar <- NULL

  oos <- matrix(Inf, nrow=nfold, ncol=argl$nlambda,
                dimnames=list(levels(foldid),names(lambda)))

  if(verb) cat("fold ")

  ## define the folddev function
  folddev <- function(k){
    require(gamlr)
    train <- which(foldid!=k)
    if(!is.null(argl$shift)) argl$shift <- shift[train]
    suppressWarnings(fit <- do.call(gamlr, 
      c(list(x=x[train,],y=y[train]), argl)))
    eta <- predict(fit, x[-train,], select=0)
    if(!is.null(argl$shift)) eta <- eta + shift[-train]

    dev <- apply(eta,2, 
      function(e) 
        mean(switch(fam, 
          "gaussian" = (e-y[-train])^2, 
          "binomial" = -2*(y[-train]*e - log(1+exp(e))),
          "poisson" = -2*(y[-train]*e - exp(e)))))

    if(fam=="poisson"){
      satnllhd <- mean(ifelse(y[-train]>0,
                      y[-train]*log(y[-train]),
                      0.0) - y[-train]) 
      dev <- dev + 2*satnllhd }
    if(verb) cat(sprintf("%s,",k))
    if(length(dev) < argl$nlambda) 
      dev <- c(dev,rep(Inf,argl$nlambda-length(dev)))
    return(dev)
  }

  # apply the folddev function
  if(!is.null(cl)){
    if (requireNamespace("parallel", quietly = TRUE)) {
      parallel::clusterExport(cl,
        c("x","y","foldid","argl","fam","verb"), 
        envir=environment())
      oos <- t(parallel::parSapply(cl,1:nfold,folddev))
    } else {
      warning("cl is not NULL, but parallle package unavailable.")
      cl <- NULL
    }
  }
  if(is.null(cl)) oos <- t(sapply(1:nfold,folddev))

  cvm <- apply(oos,2,mean)
  cvs <- apply(oos,2,sd)/sqrt(nfold-1)

  seg.min <- which.min(cvm)
  lambda.min = lambda[seg.min]

  cv1se <- (cvm[seg.min]+cvs[seg.min])-cvm
  seg.1se <- min((1:length(cvm))[cv1se>=0])
  lambda.1se = lambda[seg.1se]

  if(verb) cat("done.\n")
  out <- list(gamlr=full,
          family=fam,
          nfold=nfold,
          foldid=foldid,
          cvm=cvm,
          cvs=cvs,
          seg.min=seg.min,
          seg.1se=seg.1se,
          lambda.min=lambda.min,
          lambda.1se=lambda.1se)
  class(out) <- "cv.gamlr"
  invisible(out)
}

## S3 method functions

plot.cv.gamlr <- function(x, select=TRUE, ...){

  argl = list(...)

  argl$x <- log(x$gamlr$lambda)
  argl$y <- x$cvm
  argl$type <- "n"

  if(is.null(argl$xlab)) argl$xlab="log lambda"
  if(is.null(argl$ylab)){
    if(x$family=="gaussian") argl$ylab="mean squared error"
    else argl$ylab=sprintf("%s deviance",x$family) }
  if(is.null(argl$pch)) argl$pch=20
  if(is.null(argl$col)) argl$col=4

  cvlo <- x$cvm-x$cvs
  cvhi <- x$cvm+x$cvs

  if(is.null(argl$ylim)) 
    argl$ylim=range(c(cvlo,cvhi),finite=TRUE)
  if(is.null(argl$xlim))
    argl$xlim=range(argl$x[is.finite(argl$y)])

  suppressWarnings(do.call(plot, argl))
  segments(x0=argl$x, y0=cvlo, y1=cvhi, col="grey70")
  argl$type <- NULL
  suppressWarnings(do.call(points, argl))

  if(select){
    abline(v=log(x$lambda.min), lty=3, col="grey20")
    abline(v=log(x$lambda.1se), lty=3, col="grey20") }

  dfi <- unique(round(
    seq(1,length(argl$x),length=ceiling(length(axTicks(1))))))
  axis(3,at=argl$x[dfi], 
    labels=round(x$gamlr$df[dfi],1),tick=FALSE, line=-.5)

}

coef.cv.gamlr <- function(object, 
                          select=c("1se","min"), ...){
  seg = paste("seg",match.arg(select),sep=".")
  coef(object$gamlr, select=object[[seg]])
}

predict.cv.gamlr <- function(object, newdata,
                          select=c("1se","min"), ...){
  seg = paste("seg",match.arg(select),sep=".")
  predict.gamlr(object$gamlr, newdata, select=object[[seg]], ...)
}

summary.cv.gamlr <- function(object, ...){
  print(object)

  return(data.frame(
    lambda=object$gamlr$lambda,
    par=diff(object$gamlr$b@p)+1,
    oos.r2=1-object$cvm/object$cvm[1]))
}

print.cv.gamlr <- function(x, ...){
  cat("\n")
  cat(sprintf(
    "%d-fold %s cv.gamlr object", 
    x$nfold, x$gamlr$family))
  cat("\n\n")
}














