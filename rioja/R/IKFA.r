IKFA <- function(y, x, nFact = 5, IsPoly = FALSE, IsRot = TRUE, ccoef = 1:nFact, check.data=TRUE, lean=FALSE, ...)
{
  if (check.data) {
    if (any(apply(y, 1, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following rows:", paste(which(apply(y, 1, sum) < 1.0E-8), collapse=",")))
    if (any(apply(y, 2, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following columns:", paste(which(apply(y, 2, sum) < 1.0E-8), collapse=",")))
  }
  nFact <- max(nFact, max(ccoef))
  fit <- IKFA.fit(y=y, x=x, nFact=nFact, IsPoly=IsPoly, IsRot=IsRot, ccoef=ccoef)
  xHat <- predict.internal.IKFA(object=fit, y=y, lean=lean, ...) 
  call.fit <- as.call(list(quote(IKFA.fit), y=quote(y), x=quote(x), nFact=nFact, IsPoly=IsPoly, ccoef=ccoef, lean=FALSE))
  call.print <- match.call()
  result <- c(fit, list(fitted.values=xHat, call.fit=call.fit, call.print=call.print, x=x))
  result$cv.summary <- list(cv.method="none")
	if (!lean) 
	   result$y <- y
  class(result) <- "IKFA" 
  result
}

IKFA.fit <- function(y, x, nFact = 5, IsPoly = FALSE, IsRot = TRUE, ccoef = 1:nFact, lean=FALSE)
{ 
  nFact <- max(nFact, max(ccoef))
  y <- as.matrix(y)
	W <- y/sqrt(apply(y * y, 1, FUN = sum))
	e <- eigen(crossprod(W, W), symmetric=TRUE)
	V <- as.matrix(e$vectors[, 1:nFact])
  rownames(V) <- colnames(y)
	U <- W %*% V
#  VRot <- rotate (V, rotation="varimax", normalize=F)$rmat
  VRot=NULL
	if (IsRot) {
     VRot <- varimax (V, normalize=FALSE)$loadings
     VRot <- VRot[, 1:nFact]
     rownames(VRot) <- colnames(y)
     URot <- W %*% VRot
	   U2 <- URot
  } else {
     URot=NULL
     U2 <- U
  }
  U2 <- U2[, ccoef, drop=FALSE]
  regr <- vector("list", ncol(U2))
	for (i in 1:ncol(U2)){
	 if (IsPoly) {
    	rgr <- lm (x ~ poly(U2[, 1:i, drop=FALSE], degree=2, raw=TRUE))
    }
    else {
      rgr <- lm (x ~ U2[, 1:i, drop=FALSE])
    }
    regr[[i]] <- rgr
  }
	comm <- apply((U)^2, 1, sum)
	result <- list("score" = U, "V" = V, "score.rot"=URot, "VRot"=VRot, "comm" = comm, eig=e, "regr" = regr, "Rotated" = IsRot, "IsPoly" = IsPoly, nFact=nFact, ccoef=ccoef)
	result
}

communality <- function(object, y)
{
  y <- as.matrix(y)
	n <- nrow(y)
	W <- y/sqrt(apply(y * y, 1, FUN = sum))
  U <- W %*% object$V
	if (object$Rotated) {
    	U <- W %*% object$VRot
  }
  apply((U)^2, 1, sum)
}

predict.internal.IKFA <- function(object, y, lean=FALSE, ...)
{
  y <- as.matrix(y)
	n <- nrow(y)
	W <- y/sqrt(apply(y * y, 1, FUN = sum))
	if (object$Rotated) {
    	U <- W %*% object$VRot
  } else {
    U <- W %*% object$V
  }
  nFact <- max(object$ccoef)
  xHat <- matrix(NA, nrow=nrow(y), ncol=length(object$ccoef))
	U2 <- U[, object$ccoef, drop=FALSE]
  for (i in 1:length(object$ccoef)) {
    if (object$IsPoly) {
      # workaround for bug in poly when called with matrix with one row
      if (nrow(U2)==1) {
        U2 <- rbind(U2, U2)
        tmp <- cbind(rep(1, n), poly(U2[, 1:i, drop=FALSE], degree=2, raw=TRUE)) %*% object$regr[[i]]$coefficients
        xHat[, i] <- tmp[1, ]
      } else {
        xHat[, i] <- cbind(rep(1, n), poly(U2[, 1:i, drop=FALSE], degree=2, raw=TRUE)) %*% object$regr[[i]]$coefficients
      }
    }
    else {
      xHat[, i] <- cbind(rep(1, n), U2[, 1:i, drop=FALSE]) %*% object$regr[[i]]$coefficients
    }
  }
  rownames(xHat) <- rownames(y)
  colnames(xHat) <- paste("Fact", sprintf("%02d", object$ccoef), sep="")
  xHat
}

crossval.IKFA <- function(object, cv.method="loo", verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, ...) {
  .crossval(object=object, cv.method=cv.method, verbose=verbose, ngroups=ngroups, nboot=nboot, h.cutoff=h.cutoff, h.dist=h.dist, ...)
}

predict.IKFA <- function(object, newdata=NULL, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, ...) {
  .predict(object=object, newdata=newdata, sse=sse, nboot=nboot, match.data=match.data, verbose=verbose, ...)
}

print.IKFA <- function(x, ...) 
{
  cat("\n")
  cat("Method : Imbrie & Kipp Factor Analysis\n")
  cat("Call   : ")
  cat(paste(deparse(x$call.print), "\n\n"))
  cat(paste("No. samples        :", length(x$x), "\n"))
  cat(paste("No. species        :", nrow(x$V), "\n"))
  .print.crossval(x)
  cat("\nPerformance:\n")
  .print.performance(x)
  cat("\n")
}

performance.IKFA <- function(object, ...) {
  .performance(object, ...)
}

summary.IKFA <- function(object, full=FALSE, ...) 
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

plot.IKFA <- function(x, resid=FALSE, xval=FALSE, nFact=max(x$ccoef), xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE, add.smooth=FALSE, ...) {
  if (xval & x$cv.summary$cv.method=="none")
     stop("IKFA model does not have cross validation estimates")
  xx <- x$x
  if (resid) {
     if (xval) {
       yy <- x$predicted[, nFact] - x$x
     } else {
       yy <- residuals(x)[, nFact]
     }
  } else {
     if (xval) {
        yy <- x$predicted[, nFact]
      }  else {
       yy <- x$fitted.values[, nFact]
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

fitted.IKFA <- function(object, ...) {
  object$fitted.values
}

residuals.IKFA <- function(object, cv=FALSE, ...) {
  if (cv == FALSE)
     return (object$x - object$fitted.values)
  else {
     if (object$cv.summary$cv.method == "none")
        stop("Object does not contain cross validation results")
     return (object$residuals.cv)
  }
}

coef.IKFA <- function(object, ...) {
  object$V
}

screeplot.IKFA <- function(x, rand.test=TRUE, ...) {
  if (x$cv.summary$cv.method=="none")
     stop("IKFA model does not have cross validation estimates")
  summ <- performance(x)
  yR <- range(summ$object[, "RMSE"], summ$crossval[, "RMSE"])
  oldmar <- par("mar")
  if (rand.test) {
    mm <- oldmar
    mm[4] <- mm[4] * 2
    par(mar=mm)
  }
  plot(1:length(x$ccoef), summ$object[, "RMSE"], type="b", ylim=yR, col="black", axes=FALSE, xlab="Number of factors", ylab="RMSE")
  axis(2, las=1)
  axis(1, at=1:length(x$ccoef), labels=1:length(x$ccoef))
  box()
  lines(1:length(x$ccoef), summ$crossval[, "RMSE"], type="b", col="red")
  args <- names(as.list(match.call()))
  if ("main" %in% args) 
     title(...)
  else 
     title(main=substitute(x))
  ncomp <- nrow(summ$crossval)
  if (rand.test) {
    res <- rand.t.test(x, ...)
    rY <- range(pretty(res[-1, "delta.RMSE"]), na.rm=TRUE)
    us <- par("usr")
    us[3] <- rY[1]
    us[4] <- rY[2]
    par(usr=us)
    axis(4, las=1)
    lines(1:ncomp, res[, "delta.RMSE"], type="b", col="blue")
    mtext("% change in RMSE", 4, line=2)
    text(1:ncomp, res[, "delta.RMSE"], sprintf("%.3f", res[, "p"], 3), pos=2, xpd=NA, cex=0.8, col="blue") 
    legend("bottomleft", c("model RMSE", "x-val RMSE", "% change RMSE"), lty=1, col=c("black", "red", "blue"))
    par(mar=oldmar)
    invisible(res)
  }
  else  
    legend("topright", c("model", "x-val"), lty=1, col=c("black", "red"))
}

rand.t.test.IKFA <- function(object, n.perm=999, ...) 
{ 
  if (object$cv.summary$cv.method=="none")
     stop("Can only perform a randomisation t-test on a cross-validated model.")
  p <- performance(object)
  ncomp <- nrow(p$crossval)
  if (ncomp < 2)
     stop("Only one factor - nothing to test.")
  delta <- diff(p$crossval[, 1]) / p$crossval[1:(ncomp-1)] * 100
  e <- object$residuals.cv
  t.res <- vector("numeric",ncomp)
  t.res[] <- NA
  for (i in 2:ncomp) {
    d <- e[, i-1]^2 - e[, i]^2
    t.obs <- mean(d, na.rm=TRUE)
    t.sum <- 0
    for (j in 1:n.perm) {
      rnd <- sample(c(TRUE, FALSE), 100, replace=TRUE)
      d2 <- ifelse(rnd, abs(d), -abs(d))
      t <- mean(d2, na.rm=TRUE)
      if (t >= t.obs)
         t.sum <- t.sum + 1
    }
    t.res[i] <- t.sum / (n.perm+1)
  }
  result <- cbind(p$crossval, delta.RMSE=c(NA, delta), p=t.res)
  result
}

