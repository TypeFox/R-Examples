WA <- function(y, x, mono=FALSE, tolDW = FALSE, use.N2=TRUE, tol.cut=.01, check.data=TRUE, lean=FALSE)
{
#  if (mono) {
#     if (!require(mgcv, quietly=TRUE)) {
#        stop("Function WA with monotonic deshrinking requires package mgcv")
#     }
#  }
  y <- as.matrix(y)
  if (!is.null(dim(x)))
     x <- x[, 1]
  if (check.data) {
    if (any(apply(y, 1, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following rows:", paste(which(apply(y, 1, sum) < 1.0E-8), collapse=",")))
    if (any(apply(y, 2, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following columns:", paste(which(apply(y, 2, sum) < 1.0E-8), collapse=",")))
  }
  model <- WA.fit(y, x, mono=mono, tolDW, use.N2=use.N2, tol.cut=tol.cut, lean=lean)
  xHat <- predict.internal.WA(model, y, lean=lean)
  call.fit <- as.call(list(quote(WA.fit), y=quote(y), x=quote(x), mono=mono, tolDW=tolDW, use.N2=use.N2, tol.cut=tol.cut, lean=TRUE))
  call.print <- match.call()
	result <- c(model, list(fitted.values=xHat, call.fit=call.fit, call.print=call.print, x=x))
	if (!lean) 
	   result$y <- y
	class(result) <- c("WA")
  result$cv.summary <- list(cv.method="none")
	result
}

WA.fit <- function(y, x, mono=FALSE, tolDW=FALSE, use.N2=TRUE, tol.cut=0.01, lean=FALSE) {
  if (is.data.frame(y))
     y <- as.matrix(y)
  cS <- colSums(y, na.rm=TRUE)
  rS <- rowSums(y, na.rm=TRUE)
  u <- colSums(y * x) / cS
  bad <- is.nan(u)
  u[bad] <- NA
  u2 <- u
#  u2[bad] <- 1E8
  u2[bad] <- 0 # need to check this
  xHat <- (y %*% u2 / rS)[, 1]
  xHat[is.nan(xHat)] <- NA
  x[rS<1.0E-8] <- NA
  x.mean <- mean(x, na.rm=TRUE)
  xHat.mean <- mean(xHat, na.rm=TRUE)
	x.scaled <- x - x.mean
	xHat.scaled <- xHat - xHat.mean
  b1.inv = sum(x.scaled * xHat.scaled, na.rm=TRUE) / sum(xHat.scaled^2, na.rm=TRUE)
  b0.inv = x.mean - b1.inv*xHat.mean
  b1.cla = sum(x.scaled * xHat.scaled, na.rm=TRUE) / sum(x.scaled^2, na.rm=TRUE)
  b0.cla = xHat.mean - b1.cla * x.mean
  deshrink.coefficients <- matrix(c(b0.inv,b1.inv,b0.cla,b1.cla), ncol=2)
  if (!lean) {
    colnames(deshrink.coefficients) <- c("Inverse d/s", "Classical d/s")
    rownames(deshrink.coefficients) <- c("wa.b0", "wa.b1")
  }
  if (mono) {
    mono.mod <- .mono.fit(xHat, x)
  }
  if (tolDW) {
    tol <- vector(mode="numeric", length=ncol(y))
    tol <- colSums(y * outer(x, u, "-")^2) / cS
    tol[is.nan(tol)] <- NA
    if (use.N2) {
       N2 <- Hill.N2(y)
       tol <- tol / (1 - 1/N2)
    }
    tol[is.infinite(tol)] <- NA
    tol <- sqrt(tol)
    tol.mean <- mean(tol, na.rm=TRUE)
#    tol[tol < tol.mean/10] <- tol.mean
    tol[is.na(tol) | tol < tol.cut] <- tol.mean
    tol[is.na(tol) | tol < tol.cut] <- 1
    tol2 <- tol^2
    rS.t <- rowSums(sweep(y, 2, tol2, "/"), na.rm=TRUE)
    tol[bad] <- NA
    u.tol <- u2/tol2
#    u.tol[bad] <- 1E8
    xHat.t <- (y %*% u.tol / rS.t)[, 1]
    xHat.t[is.nan(xHat.t)] <- NA
    xHat.t.mean <- mean(xHat.t, na.rm=TRUE)
	  xHat.t.scaled <- xHat.t - xHat.t.mean

    b1.inv.t = sum(x.scaled * xHat.t.scaled, na.rm=TRUE) / sum(xHat.t.scaled^2, na.rm=TRUE)
    b0.inv.t = x.mean - b1.inv.t*xHat.t.mean
    b1.cla.t = sum(x.scaled * xHat.t.scaled, na.rm=TRUE) / sum(x.scaled^2, na.rm=TRUE)
    b0.cla.t = xHat.t.mean - b1.cla.t * x.mean
    coef <- cbind(u, tol)
    colnames(coef) <- c("Optima", "Tolerances")
    deshrink.coefficients.t <- matrix(c(b0.inv.t,b1.inv.t,b0.cla.t,b1.cla.t), ncol=2)
    rownames(deshrink.coefficients.t) <- c("wa.tol.b0", "wa.tol.b1")
    deshrink.coefficients <- rbind(deshrink.coefficients, deshrink.coefficients.t)
    if (mono) {
      mono.mod.td <- .mono.fit(xHat.t, x)
    }
  }
  else {
    coef <- as.matrix(u)
    colnames(coef) <- "Optima"
  }
  result <- list(coefficients=coef, deshrink.coefficients=deshrink.coefficients, mono=mono, tolDW=tolDW)
  if (mono) {
     if (tolDW) {
        result$mono.mod <- list(WA=mono.mod, tolDW=mono.mod.td)
     } else {
        result$mono.mod <- list(WA=mono.mod)
     }
  }
  result
}

predict.internal.WA <- function(object, y, lean=FALSE, ...) {
  if (!lean) {
    if (nrow(object$coefficients) != ncol(y))
       stop("Number of taxa does not match that in model")
    if (any(colnames(y) != rownames(object$coefficients)))
       stop("Taxon names do not match those in model")
  }
  u <- object$coefficients[, 1]
  u.missing <- is.na(u)
  u[u.missing] <- 0
	xHat1 <- y %*% u / rowSums(y[, !u.missing, drop=FALSE])
  xHat.inv <- object$deshrink.coefficients[1, 1] + object$deshrink.coefficients[2, 1] * xHat1
  xHat.cla <- (xHat1 - object$deshrink.coefficients[1, 2]) / object$deshrink.coefficients[2, 2]
	xHat <- cbind(xHat.inv, xHat.cla)
  if (!lean)
	   colnames(xHat) <- c("WA.inv", "WA.cla")
	xHat[is.nan(xHat)] <- NA
	if (object$mono) {
     xHat.mono <- .mono.predict(object$mono.mod$WA, xHat1)
     xHat <- cbind(xHat, xHat.mono)
     colnames(xHat)[3] <- "WA.m"
  }
  if (object$tolDW) {
    tol2 <- object$coefficients[, 2]^2
    rS.t <- rowSums(sweep(y, 2, tol2, "/"), na.rm=TRUE)
    u.tol <- u/tol2
    u.tol[u.missing] <- 0
    xHat.t <- y %*% u.tol / rS.t
    xHat.t[is.nan(xHat.t)] <- NA
    xHat.inv.t <- object$deshrink.coefficients[3, 1] + object$deshrink.coefficients[4, 1] * xHat.t
    xHat.cla.t <- (xHat.t - object$deshrink.coefficients[3, 2]) / object$deshrink.coefficients[4, 2]
    xHat.t2 <- cbind(xHat.inv.t, xHat.cla.t)
    if (!lean)
    	colnames(xHat.t2) <- c("WA.inv.tol", "WA.cla.tol")
    xHat <- cbind(xHat, xHat.t2)
  	if (object$mono) {
       xHat.mono.t <- .mono.predict(object$mono.mod$tolDW, xHat.t)
       xHat <- cbind(xHat, xHat.mono.t)
       colnames(xHat)[6] <- "WA.m.tol"
    }
  }
	xHat
}

crossval.WA <- function(object, cv.method="loo", verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, ...) {
  .crossval(object=object, cv.method=cv.method, verbose=verbose, ngroups=ngroups, nboot=nboot, h.cutoff=h.cutoff, h.dist=h.dist, ...)
}

predict.WA <- function(object, newdata=NULL, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, ...) {
  .predict(object=object, newdata=newdata, sse=sse, nboot=nboot, match.data=TRUE, verbose=verbose, ...)
}

print.WA <- function(x, ...) 
{
  cat("\n")
  cat("Method : Weighted Averaging\n")
  cat("Call   : ")
  cat(paste(deparse(x$call.print), "\n\n"))
  cat(paste("Tolerance DW       :", ifelse(x$tolDW, "Yes", "No"), "\n"))
  cat(paste("Monotonic deshrink :", ifelse(x$mono, "Yes", "No"), "\n"))
  cat(paste("No. samples        :", length(x$x), "\n"))
  cat(paste("No. species        :", nrow(x$coefficients), "\n"))
  .print.crossval(x)
  cat("Deshrinking regression coefficients:\n")
  print(round(x$deshrink.coefficients, 4))
  cat("\nPerformance:\n")
  .print.performance(x)
  cat("\n")
}

performance.WA <- function(object, ...) {
  .performance(object, ...)
}

summary.WA <- function(object, full=FALSE, ...) 
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

plot.WA <- function(x, resid=FALSE, xval=FALSE, tolDW=FALSE, deshrink="inverse", xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE, add.smooth=FALSE, ...) {
  DESHRINK <- c("inverse", "classical", "mono")
  deshrink <- pmatch(deshrink, DESHRINK)
 
  if (is.na(deshrink))
     stop("Unknown deshrinking method")

  if (deshrink == 3 & !x$mono) 
     stop("WA model does not have monotonic spline deshrinking estimates")
  
  if (tolDW) {
    if (!x$tolDW) {
       stop("WA model does not have tolerance downweighting estimates")
     } else {
       if (x$mono) {
         deshrink <- deshrink + 3
       } else {
         deshrink <- deshrink + 2
       }
     }
  }
  if (xval & x$cv.summary$cv.method=="none")
     stop("WA model does not have cross validation estimates")
  xx <- x$x
  if (resid) {
     if (xval) {
       yy <- x$predicted[, deshrink] - x$x
     } else {
       yy <- residuals(x)[, deshrink]
     }
  } else {
     if (xval) {
        yy <- x$predicted[, deshrink]
      }  else {
        yy <- x$fitted.values[, deshrink]
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

fitted.WA <- function(object, ...) {
  object$fitted.values
}

residuals.WA <- function(object, cv=FALSE, ...) {
  if (cv == FALSE) {
     return(object$x - object$fitted.values)
  } else {
     if (object$cv.summary$cv.method == "none")
        stop("Object does not contain cross validation results")
     return(object$residuals.cv)
  }
}

coef.WA <- function(object, ...) {
  object$coefficients
}

rand.t.test.WA <- function(object, n.perm=999, ...) 
{ 
  if (object$cv.summary$cv.method=="none")
     stop("Can only perform a randomisation t-test on a cross-validated model.")
  p <- performance(object)
  if (!object$tolDW & !object$mono)
     stop("Model does not contain tolerance downweighting or monotonic deshrinking - nothing to test.")
  mono <- ifelse(object$mono, 1, 0)
  tolDW <- ifelse(object$tolWD, 2, 0)
  e <- object$residuals.cv
  delta1 <- delta2 <- d1 <- d2 <- t1.obs <- t2.obs <- t1.sum <- t2.sum <- 0
  if (object$tolDW) {
     delta1 <- (p$crossval[3+mono, 1] - p$crossval[1, 1]) / p$crossval[1, 1] * 100
     delta2 <- (p$crossval[4+mono, 1] - p$crossval[2, 1]) / p$crossval[2, 1] * 100
     d1 <- e[, 1]^2 - e[, 3+mono]^2
     d2 <- e[, 2]^2 - e[, 4+mono]^2
     t1.obs <- mean(d1, na.rm=TRUE)
     t2.obs <- mean(d2, na.rm=TRUE)
     t1.sum <- 0
     t2.sum <- 0
     if (object$mono) {
        delta4 <- (p$crossval[6, 1] - p$crossval[4, 1]) / p$crossval[4, 1] * 100
        d4 <- e[, 4]^2 - e[, 6]^2
        t4.obs <- mean(d4, na.rm=TRUE)
        t4.sum <- 0
     }
  }
  if (object$mono) {
     delta3 <- (p$crossval[3, 1] - p$crossval[1, 1]) / p$crossval[1, 1] * 100
     d3 <- e[, 1]^2 - e[, 3]^2
     t3.obs <- mean(d3, na.rm=TRUE)
     t3.sum <- 0
  }
  for (j in 1:n.perm) {
     rnd <- sample(c(TRUE, FALSE), 100, replace=TRUE)
     if (object$tolDW) {
       d12 <- ifelse(rnd, abs(d1), -abs(d1))
       d22 <- ifelse(rnd, abs(d2), -abs(d2))
       t1 <- mean(d12, na.rm=TRUE)
       t2 <- mean(d22, na.rm=TRUE)
       if (t1 >= t1.obs)
         t1.sum <- t1.sum + 1
       if (t2 >= t2.obs)
         t2.sum <- t2.sum + 1
       if (object$mono) {
         d14 <- ifelse(rnd, abs(d4), -abs(d4))
         t4 <- mean(d14, na.rm=TRUE)
         if (t4 >= t4.obs)
           t4.sum <- t4.sum + 1
      }
     }
     if (object$mono) {
       d13 <- ifelse(rnd, abs(d3), -abs(d3))
       t3 <- mean(d13, na.rm=TRUE)
       if (t3 >= t3.obs)
         t3.sum <- t3.sum + 1
     }
  }
  
  if (object$tolDW) {
    t1.res <- t1.sum / (n.perm+1)
    t2.res <- t2.sum / (n.perm+1)
    if (object$mono) {
      t3.res <- t3.sum / (n.perm+1)
      t4.res <- t4.sum / (n.perm+1)
      result <- cbind(p$crossval, delta.RMSE=c(NA, NA, delta3, delta1, delta2, delta4), p=c(NA, NA, t3.res, t1.res, t2.res, t4.res))
    } else {
      result <- cbind(p$crossval, delta.RMSE=c(NA, NA, delta1, delta2), p=c(NA, NA, t1.res, t2.res))
    }
  } else {
    t3.res <- t3.sum / (n.perm+1)
    result <- cbind(p$crossval, delta.RMSE=c(NA, NA, delta3), p=c(NA, NA, t3.res))
  }
  result
}

.mono.fit<-function(x,y) {
   dat <- data.frame(x=x,y=y)
   f.ug <- gam(y~s(x,k=10,bs="cr"));
  # Create Design matrix, constraints etc. for monotonic spline....
   sm<-smoothCon(s(x,k=10,bs="cr"),dat,knots=NULL)[[1]]
   Fm<-mono.con(sm$xp);   # get constraints
   G<-list(X=sm$X,C=matrix(0,0,0),sp=f.ug$sp,p=sm$xp,y=y,w=y*0+1)
   G$Ain<-Fm$A;G$bin<-Fm$b;G$S<-sm$S;G$off<-0
   p<-pcls(G);  # fit spline (using s.p. from unconstrained fit)
   list(p=p,sm=sm)
}

.mono.predict<-function(mod, newdata) {
   Predict.matrix(mod$sm,data.frame(x=newdata))%*%mod$p
}
