MLRC <- function(y, x, check.data=TRUE, lean=FALSE, n.cut=5, verbose=TRUE, ...)
{
  if (check.data) {
    if (any(apply(y, 1, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following rows:", paste(which(apply(y, 1, sum) < 1.0E-8), collapse=",")))
    if (any(apply(y, 2, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following columns:", paste(which(apply(y, 2, sum) < 1.0E-8), collapse=",")))
    if(n.cut < 5 & any(apply(y>0, 2, sum) < 5))
       warning("Trying to fit responses to some taxa with less than 5 occurrences - results may be unreliable")
  }
  if (any(y>1) | any (y<0))
     stop("Species data must be proportions between 0 and 1")
  fit <- MLRC.fit(y=y, x=x, lean=lean, n.cut=n.cut, verbose=verbose, ...)
  xHat <- predict.internal.MLRC(object=fit, y=y, lean=lean, ...)
  call.print <- match.call()
  call.fit <- as.call(list(quote(MLRC.fit), y=quote(y), x=quote(x), lean=FALSE))
  result <- c(fit, list(fitted.values=xHat, call.fit=call.fit, call.print=call.print, x=x))
  result$cv.summary <- list(cv.method="none")
	if (!lean) 
	   result$y <- y
  class(result) <- "MLRC" 
  result
}

MLRC.fit <- function(y, x, n.cut=2, use.glm = FALSE, max.iter=50, lean=FALSE, verbose=TRUE, ...)
{ 
  glr <- function(x, e) {
    gfit <- glm.fit(e, x, family = quasibinomial(link=logit), ...)
    if (gfit$converged)
       return(gfit$coefficients)
    else
       return(c(NA, NA, NA))
  }
  skip <- colSums(y > 0) < n.cut
  if (use.glm) {
#    glr <- function(x, e) {
#	 	   gfit <- glm(x ~ e + I(e^2), family = quasibinomial(link=logit), ...)
#	 	   if (gfit$converged)
#		      return(gfit$coefficients)
#      else
#          return(c(NA, NA, NA))
#    }
    lp <- cbind(rep(1, nrow(y)), x, x^2)
	  beta <- apply(y[, !skip], 2, glr, e=lp)
    BETA <- matrix(NA, nrow = 3, ncol = ncol(y))
    BETA[, !skip] <- beta
    rownames(beta) <- c("b0", "b1", "b2")
	  return (list(coefficients=t(BETA), meanX=mean(x, na.rm=TRUE)))
  } else {
    res <- .Call("MLRC_regress", as.matrix(y[, !skip]), as.matrix(x), as.integer(max.iter), as.integer(verbose), NAOK=TRUE, PACKAGE="rioja")  
    beta <- matrix(res$Beta, ncol=3)
    BETA <- matrix(NA, ncol = 3, nrow = ncol(y))
    BETA[!skip, ] <- beta
    IBETA <- vector("integer", length=ncol(y))
    IBETA[] <- NA
    IBETA[!skip] <- res$IBeta
    rownames(BETA) <- colnames(y)
    colnames(BETA) <- c("b0", "b1", "b2")
    list(coefficients=BETA, meanX=mean(x, na.rm=TRUE), IBeta=IBETA, n.cut=n.cut)
  }
}


predict.internal.MLRC <- function(object, y, lean=FALSE, ...)
{
	coef <- object$coefficients
	if (!lean) {
	   if (nrow(object$coefficients) != ncol(y))
	      stop("Number of columns different in y, beta in predict.internal.MLRC")
	}
  xHat <- .Call("MLRC_predict", as.matrix(y), as.matrix(object$coefficients), as.double(object$meanX), NAOK=TRUE, PACKAGE="rioja")
  xHat <- as.matrix(xHat, ncol=1)
  colnames(xHat) <- "MLRC"
  rownames(xHat) <- rownames(y)
  xHat
}

crossval.MLRC <- function(object, cv.method="loo", verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, ...) {
  .crossval(object=object, cv.method=cv.method, verbose=verbose, ngroups=ngroups, nboot=nboot, h.cutoff=h.cutoff, h.dist=h.dist, ...)
}

predict.MLRC <- function(object, newdata=NULL, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, ...) {
   if (!is.null(newdata))
      if (any(newdata < 0) | any(newdata > 1))
         stop("newdata must be proportions between 0 and 1")
  .predict(object=object, newdata=newdata, sse=sse, nboot=nboot, match.data=match.data, verbose=verbose, ...)
}

performance.MLRC <- function(object, ...) {
  .performance(object, ...)
}

print.MLRC <- function(x, ...) 
{
  cat("\n")
  cat("Method : Maximum Likelihood using Response Curves \n")
  cat("Call   : ")
  cat(paste(deparse(x$call.print), "\n\n"))
  cat(paste("No. samples        :", length(x$x), "\n"))
  cat(paste("No. species        :", nrow(x$coefficients), "\n"))
  .print.crossval(x)
  cat("\nPerformance:\n")
  .print.performance(x)
  cat("\n")
}

summary.MLRC <- function(object, full=FALSE, ...) 
{
  print(object, ...)
  if (object$cv.summary$cv.method == "none")
    fitted <- as.data.frame(object$fitted.values)     
  else
    fitted <- as.data.frame(object$fitted.values, object$predicted)     
  cat("\nFitted values\n")
  if (full) {
     print(fitted)
     cat("\nSpecies coefficients\n")
     print(data.frame(object$coefficients))
  } else {
     print(dot(fitted))
     cat("\nSpecies coefficients\n")
     print(dot(data.frame(object$coefficients)))
  }
}

plot.MLRC <- function(x, resid=FALSE, xval=FALSE, xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE, add.smooth=FALSE, ...) {
  if (xval & x$cv.summary$cv.method=="none")
     stop("MLRC model does not have cross validation estimates")
  xx <- x$x
  if (resid) {
     if (xval) {
       yy <- x$predicted[, 1]
     } else {
       yy <- residuals(x)[, 1]
     }
  } else {
     if (xval) {
        yy <- x$predicted[, 1]
      }  else {
       yy <- x$fitted.values[, 1]
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

fitted.MLRC <- function(object, ...) {
  object$fitted.values
}

residuals.MLRC <- function(object, cv=FALSE, ...) {
  if (cv == FALSE)
     return (object$x - object$fitted.values)
  else {
     if (object$cv.summary$cv.method == "none")
        stop("Object does not contain cross validation results")
     return (object$residuals.cv)
  }
}

coef.MLRC <- function(object, ...) {
  object$coefficients
}

#predict.internal.MLRC <- function(object, y, lean=FALSE, ...)
#{
#  y <- as.matrix(y)
#	nnn <- nrow(y)
#	xresp <- object$xSearch
#	yresp <- object$resp
#	nn <- length(xresp)
#	p <- log(yresp)
#	ppp <- log(1-yresp)
# LL.res <- as.matrix(p) %*% t(y) + as.matrix(ppp) %*% t(1.0-y)
#  LL.res[is.na(LL.res)] <- -1.0E10
#	xHat <- xresp[apply(LL.res, 2, order)[nn, ]]
#  xHat <- as.matrix(xHat, ncol=1)
#  colnames(xHat) <- "MLRC"
#  rownames(xHat) <- rownames(y)
#  xHat
#}

#MLRC.fit <- function(y, x, xSearch, lean=FALSE)
#{ 
#  glr <- function(x, e, xSearch) {
#		 gfit <- glm(x ~ e + I(e^2), family = quasibinomial(link=logit))
#		 gfit$coefficients
#		 predict.glm(gfit, data.frame(e=xSearch), type="response")
#  }
#	resp <- apply(y, 2, glr, e=x, xSearch)
#	result <- list(resp=resp, xSearch=xSearch)
#}

