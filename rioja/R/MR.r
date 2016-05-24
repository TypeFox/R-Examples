MR <- function(y, x, check.data=TRUE, lean=FALSE, ...)
{
  if (check.data) {
    if (any(apply(y, 1, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following rows:", paste(which(apply(y, 1, sum) < 1.0E-8), collapse=",")))
    if (any(apply(y, 2, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following columns:", paste(which(apply(y, 2, sum) < 1.0E-8), collapse=",")))
  }
  fit <- MR.fit(y=y, x=x)
  xHat <- predict.internal.MR(object=fit, y=y, lean=lean, ...) 
  call.fit <- as.call(list(quote(MR.fit), y=quote(y), x=quote(x), lean=FALSE))
  call.print <- match.call()
  result <- list(lm=fit, fitted.values=xHat, call.fit=call.fit, call.print=call.print, x=x)
  result$cv.summary <- list(cv.method="none")
	if (!lean) 
	   result$y <- y
  class(result) <- "MR" 
  result
}

MR.fit <- function(y, x, lean=FALSE)
{ 
  rgr <- lm (x ~ ., data=y)
  rgr
}

predict.internal.MR <- function(object, y, ...)
{
  y <- as.data.frame(y)
  xHat <- predict.lm(object, y)
  xHat <- matrix(xHat, ncol=1)
  rownames(xHat) <- rownames(y)
  colnames(xHat) <- "xHat"
  xHat
}

crossval.MR <- function(object, cv.method="loo", verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, ...) {
  .crossval(object=object, cv.method=cv.method, verbose=verbose, ngroups=ngroups, nboot=nboot, h.cutoff=h.cutoff, h.dist=h.dist, ...)
}

predict.MR <- function(object, newdata=NULL, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, ...) {
  .predict(object=object, newdata=newdata, sse=sse, nboot=nboot, match.data=match.data, verbose=verbose, ...)
}

print.MR <- function(x, ...) 
{
  cat("\n")
  cat("Method : Multiple regression\n")
  cat("Call   : ")
  cat(paste(deparse(x$call.print), "\n\n"))
  cat(paste("No. samples        :", length(x$x), "\n"))
  cat(paste("No. species        :", length(x$lm$coefficients)-1, "\n"))
  .print.crossval(x)
  cat("\nPerformance:\n")
  .print.performance(x)
  cat("\n")
}

performance.MR <- function(object, ...) {
  .performance(object, ...)
}

summary.MR <- function(object, full=FALSE, ...) 
{
  print(object, ...)
  if (object$cv.summary$cv.method == "none")
    fitted <- as.data.frame(object$fitted.values)     
  else
    fitted <- as.data.frame(object$fitted.values, object$predicted)     
  cat("\nFitted values\n")
  if (full) {
     print(fitted)
  } else {
     print(dot(fitted))
  }
}

plot.MR <- function(x, resid=FALSE, xval=FALSE, xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE, add.smooth=FALSE, ...) {
  if (xval & x$cv.summary$cv.method=="none")
     stop("IKFA model does not have cross validation estimates")
  xx <- x$x
  if (resid) {
     if (xval) {
       yy <- x$predicted - x$x
     } else {
       yy <- residuals(x)
     }
  } else {
     if (xval) {
        yy <- x$predicted
      }  else {
       yy <- x$fitted.values
      }
  }
  if (missing(ylim)) {
     if (resid) {
       ylim <- range(yy)
     } else {
       ylim <- range(yy, x$x)
     }
  }
  if (missing(xlim))
     xlim <- range(xx, x$x)
  plot(xx, yy, ylim=ylim, xlim=xlim, xlab=xlab, ylab=ylab, las=1, ...)
  if (add.ref) {
     if (resid)
       abline(h=0, col="grey")
     else
       abline(0,1, col="grey")
  }
  if (add.smooth) {
     lines(lowess(xx, yy), col="red")
  }
}

fitted.MR <- function(object, ...) {
  object$fitted.values
}

residuals.MR <- function(object, cv=FALSE, ...) {
  if (cv == FALSE)
     return (object$x - object$fitted.values)
  else {
     if (object$cv.summary$cv.method == "none")
        stop("Object does not contain cross validation results")
     return (object$residuals.cv)
  }
}

coef.MR <- function(object, ...) {
  object$lm$coefficients
}

