MAT <- function(y, x, dist.method="sq.chord", k=5, lean=TRUE)
{
  call.fit <- as.call(list(quote(MAT), y=quote(y), x=quote(x), dist.method=dist.method, k=k, lean=lean))
  x1 <- as.numeric(x)
  n <- 2
  diss <- as.matrix(paldist(y, dist.method=dist.method))
  ind <- apply(diss, 2, order)
  dist.n <- t(apply(diss, 2, sort)[n:(n+k-1), , drop=FALSE])
  rownames(dist.n) <- rownames(y)
  colnames(dist.n) <- paste("N", sprintf("%02d", 1:k), sep="")
  x.n <- t(matrix(x1[ind[n:(n+k-1), , drop=FALSE]], nrow=k))
  rownames(x.n) <- rownames(y)
  colnames(x.n) <- paste("N", sprintf("%02d", 1:k), sep="")
  xHat <- matrix(NA, nrow=nrow(y), ncol=k*2)
  for (i in 1:k) {
    xHat[, i] <- apply(x.n[ ,1:i, drop=FALSE], 1, mean, na.rm=TRUE)
    xHat[, i+k] <- rowSums((x.n / dist.n)[ ,1:i, drop=FALSE], na.rm=TRUE) / rowSums(1/dist.n[ ,1:i, drop=FALSE])
  }
#  xHat <- cbind(xHat=xHat, xHat.wm=xHat.wm)
  rownames(xHat) <- rownames(y)
  colnames(xHat) <- c(paste("N", sprintf("%02d", 1:k), sep=""), paste("N", sprintf("%02d", 1:k), ".wm", sep=""))
  sd <- apply(x.n, 1, sd, na.rm=TRUE)
  minD <- apply(dist.n, 1, min, na.rm=TRUE)
  diagn  <- data.frame(Stdev=sd, minD=minD)  
  nms <- t(matrix(rownames(y)[ind[n:(n+k-1), , drop=FALSE]], nrow=k))
  rownames(nms) <- rownames(y)
  colnames(nms) <- paste("N", sprintf("%02d", 1:k), sep="")
  rownames(diss) <- rownames(y)
  colnames(diss) <- rownames(y)
  call.print <- match.call()
  result <- list(call=call.fit, call.print=call.print, fitted.values=xHat, diagnostics=diagn, dist.n=dist.n, x.n=x.n, match.name=nms, x=x1, dist.method=dist.method, k=k, y=y)
  if (!lean)
     result <- c(result, list(dist=diss))
  result$cv.summary <- list(cv.method="none")
	class(result) <- "MAT"
	result
}

predict.MAT <- function(object, newdata=NULL, k=object$k, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, lean=TRUE, ...) 
{
  if (k < 1 | k > object$k)
    stop("k out of range")
  if (is.null(newdata)) {
     return(object$fitted.values[, c(k, k+object$k), drop=FALSE])
  }
  if (match.data) {
    d <- Merge(object$y, newdata, split=TRUE)
    y1 <- as.matrix(d[[1]])
    y2 <- as.matrix(d[[2]])
  } else {
    if (ncol(object$y) != ncol(newdata)) 
       stop("Number of taxa does not match between datasets")
    if (any(colnames(object$y) != colnames(newdata)))
       stop("Taxon names do not match between datasets")
    y1 <- object$y
    y2 <- newdata
  }
  if (nrow(y1) != nrow(y2)) {
     n = 1
  } else if (any(rownames(y1) == rownames(y2))) {
     n = 2
  } else {
     n = 1
  }
  x1 <- object$x
#  k <- object$k
  diss <- paldist2(y1, y2, dist.method=object$dist.method)
  ind <- apply(diss, 2, order)
  dist.n <- t(apply(diss, 2, sort)[n:(n+k-1), , drop=FALSE])
  rownames(dist.n) <- rownames(y2)
  colnames(dist.n) <- paste("N", sprintf("%02d", 1:k), sep="")
  x.n <- t(matrix(x1[ind[n:(n+k-1), , drop=FALSE]], nrow=k))
  rownames(x.n) <- rownames(y2)
  colnames(x.n) <- paste("N", sprintf("%02d", 1:k), sep="")
  xHat <- apply(x.n, 1, mean, na.rm=TRUE)
  xHat.wm <- rowSums(x.n / dist.n, na.rm=TRUE) / rowSums(1/dist.n)
  xHat <- cbind(MAT=xHat, MAT.wm=xHat.wm)
  rownames(xHat) <- rownames(y2)
  sd <- apply(x.n, 1, sd, na.rm=TRUE)
  minD <- apply(dist.n, 1, min, na.rm=TRUE)
  diagn  <- data.frame(Stdev=sd, minD=minD)  
  nms <- t(matrix(rownames(y1)[ind[n:(n+k-1), , drop=FALSE]], nrow=k))
  rownames(nms) <- rownames(y2)
  colnames(nms) <- paste("N", sprintf("%02d", 1:k), sep="")
  rownames(diss) <- rownames(y1)
  colnames(diss) <- rownames(y2)
  result <- list(k=k, fit=xHat, diagnostics=diagn, dist.n=dist.n, x.n=x.n, match.name=nms)
  if (!lean)
     result <- c(result, list(dist=diss))
  if (sse) {
   feedback <- ifelse(is.logical(verbose), 50, as.integer(verbose))
   if (is.null(object$dist))
      stop("No distances: refit original model using \"lean=FALSE\"")
    nsam <- nrow(object$y)
    nsam.new <- nrow(newdata)
    nest <- 2
    res2 <- array(dim=c(nsam, nest, nboot))
    res2.new <- array(dim=c(nsam.new, nest, nboot))
#    .set.rand.seed(100)
    for (i in 1:nboot) {
#      o <- apply(data.frame(rep(nsam, nsam)), 1, .get.rand) + 1
      o <- sample(nsam, replace=TRUE)
      out <- (1:nsam)[-unique(o)]
      x <- object$x[o]

      diss.m <- object$dist[o, out]
      ind <- apply(diss.m, 2, order)
      dist.n <- t(apply(diss.m, 2, sort)[1:k, , drop=FALSE])
      x.n <- t(matrix(x[ind[1:k, , drop=FALSE]], nrow=k))
      res2[out, 1, i] <- apply(x.n, 1, mean, na.rm=TRUE)
      res2[out, 2, i] <- rowSums((x.n / dist.n), na.rm=TRUE) / rowSums(1/dist.n)
      
      diss.f <- diss[o, ]
      ind <- apply(diss.f, 2, order)
      dist.n <- apply(diss.f, 2, sort)[1:k, , drop=FALSE]
      x.n <- matrix(x[ind[1:k, , drop=FALSE]], nrow=k)
      res2.new[, 1, i] <- apply(x.n, 2, mean, na.rm=TRUE)
      res2.new[, 2, i] <- colSums((x.n / dist.n), na.rm=TRUE) / colSums(1/dist.n)
      if (verbose) {
          if (i %% feedback == 0) {
            cat (paste("Bootstrap sample", i, "\n"))
            flush.console()
          }
      }
    }      
    xHat.boot <- apply(res2, c(1,2), mean, na.rm=TRUE)
    xHat.new.boot <- apply(res2.new, c(1,2), mean, na.rm=TRUE)
    colnames(xHat.new.boot) <- colnames(result$fit)
    rownames(xHat.new.boot) <- rownames(newdata)
    v1.boot <- apply(res2.new, c(1,2), sd, na.rm=TRUE)
    v2.boot <- apply(object$x-xHat.boot, 2, .rmse)
    colnames(v1.boot) <- colnames(result$fit)
    rownames(v1.boot) <- rownames(newdata)
    SEP.boot <- sqrt(sweep(v1.boot^2, 2, v2.boot^2, "+"))
    colnames(SEP.boot) <- colnames(result$fit)
    rownames(SEP.boot) <- rownames(newdata)
    result <- c(result, list(fit.boot=xHat.new.boot, v1.boot=v1.boot, v2.boot=v2.boot, SEP.boot=SEP.boot))
  }   
  result
}

crossval.MAT <- function(object, k=object$k, cv.method="lgo", verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, ...)
{
  if (k < 1 | k > object$k)
    stop("k out of range")
  METHODS <- c("lgo", "bootstrap")
  cv.method <- pmatch(cv.method, METHODS)
  if (is.na(cv.method))
     stop("Unknown cross-validation method")
  nsam <- length(object$x)
  nres <- ncol(object$fitted.values)
  result <- matrix(nrow=nsam, ncol=nres)
  object$cv.summary$cv.method=METHODS[cv.method]
  feedback <- ifelse(is.logical(verbose), 50, as.integer(verbose))
  k <- object$k
  if (is.null(object$dist))
     stop("No distances: refit original model using \"lean=FALSE\"")
  if (cv.method == 1) {
    if (length(ngroups) > 1) {
       grps <- as.integer(factor(ngroups))
       ngroups <- length(unique(ngroups))
       o <- 1:nsam
    }
    else {
      o <- sample(nsam)
      grps <- rep(1:ngroups, length.out=nsam) 
    }
    print(ngroups)
    for (i in 1:ngroups) {
      out <- o[grps==i]
      x <- object$x[-out, drop=FALSE]
      diss <- object$dist[-out, out, drop=FALSE]
      ind <- apply(diss, 2, order)
      dist.n <- t(apply(diss, 2, sort)[1:k, , drop=FALSE])
      x.n <- t(matrix(x[ind[1:k, , drop=FALSE]], nrow=k))
      xHat <- matrix(NA, nrow=length(out), ncol=k*2)
      for (j in 1:k) {
        xHat[, j] <- apply(x.n[ ,1:j, drop=FALSE], 1, mean, na.rm=TRUE)
        xHat[, j+k] <- rowSums((x.n / dist.n)[ ,1:j, drop=FALSE], na.rm=TRUE) / rowSums(1/dist.n[ ,1:j, drop=FALSE])
      }
      result[out, ] <- xHat
      if (verbose) {
          cat (paste("Leavout group", i, "\n"))
          flush.console()
      }
    }
    object$cv.summary$ngroups=ngroups
  } else if (cv.method == 2) {
    res2 <- array(dim=c(nsam, nres, nboot))
#    .set.rand.seed(100)
    for (i in 1:nboot) {
#      o <- apply(data.frame(rep(nsam, nsam)), 1, .get.rand) + 1
      o <- sample(nsam, replace=TRUE)
      out <- (1:nsam)[-unique(o)]
      x <- object$x[o]
      diss <- object$dist[o, out]
      ind <- apply(diss, 2, order)
      dist.n <- t(apply(diss, 2, sort)[1:k, ])
      x.n <- t(matrix(x[ind[1:k, ]], nrow=k))
      for (j in 1:k) {
        res2[out, j, i] <- apply(x.n[ ,1:j, drop=FALSE], 1, mean, na.rm=TRUE)
        res2[out, j+k, i] <- rowSums((x.n / dist.n)[ ,1:j, drop=FALSE], na.rm=TRUE) / rowSums(1/dist.n[ ,1:j, drop=FALSE])
      }
      if (verbose) {
          if (i %% feedback == 0) {
            cat (paste("Bootstrap sample", i, "\n"))
            flush.console()
          }
      }
    }
    result <- apply(res2, c(1,2), mean, na.rm=TRUE)
    rownames(result) <- rownames(object$fitted.values)
    MS <- apply((object$x-res2)^2, c(1,2), mean, na.rm=TRUE)
    RMSE.boot <- sqrt(apply(MS, 2, mean, na.rm=TRUE))
    object$cv.summary$nboot=nboot
    object$cv.summary$RMSE.boot <- RMSE.boot
  } 
  colnames(result) <- colnames(object$fitted.values)
  object$predicted=result
  object$residuals.cv=result-object$x
  object
}

print.MAT <- function(x, ...) {
  cat("\n")
  cat("Method : Modern Analogue Technique\n")
  cat("Call   : ")
  cat(paste(deparse(x$call), "\n\n"))
  cat(paste("Distance :", x$dist.method, "\n"))
  cat(paste("No. samples        :", length(x$x), "\n"))
  cat(paste("No. species        :", ncol(x$y), "\n"))
  .print.crossval(x)
  cat("\nPerformance:\n")
  .print.performance(x)
  cat("\n")
}

residuals.MAT <- function(object, cv=FALSE, ...) {
  if (cv == FALSE)
     return (object$x - object$fitted.values)
  else {
     if (object$cv.summary$cv.method == "none")
        stop("Object does not contain cross validation results")
     return (object$residuals.cv)
  }
}

performance.MAT <- function(object, ...) {
  .performance(object, ...)
}

summary.MAT <- function(object, full=FALSE, ...) 
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

plot.MAT <- function(x, resid=FALSE, xval=FALSE, k=5, wMean=FALSE, xlab="", ylab="", ylim=NULL, xlim=NULL, add.ref=TRUE, add.smooth=FALSE, ...) {
  if (k < 1 | k > x$k)
     stop("k out of bounds")
  if (wMean) 
     k <- k + x$k
  if (xval & x$cv.summary$cv.method=="none")
     stop("MAT model does not have cross validation estimates")
  xx <- x$x
  if (resid) {
     if (xval) {
       yy <- x$predicted[, k] - x$x
     } else {
       yy <- residuals(x)[, k]
     }
  } else {
     if (xval) {
        yy <- x$predicted[, k]
      }  else {
        yy <- x$fitted.values[, k]
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

fitted.MAT <- function(object, ...) {
  object$fitted.values
}

screeplot.MAT <- function(x, ...) {
  summ <- performance(x)
  yR <- range(summ$object[ ,"RMSE"])
  plot(1:x$k, summ$object[1:x$k, "RMSE"], type="b", ylim=yR, col="black", ylab="RMSE", xlab="Number of analogues")
  lines(1:x$k, summ$object[(x$k+1):(2*x$k), "RMSE"], type="b", col="red")
  legend("topright", c("mean", "w-mean"), lty=1, col=c("black", "red"))
}


# 1    " Euclidean Distance          ",
# 2    " Squared Euclidean Distance  ",
# 3    " Mean Euclidean distance     ",
# 4    " Absolute Distance           ",
# 5    " Mean Absolute Distance      ",
# 6    " Percent Dissimilarity (B&C) ",
# 7    " Canberra Metric             ",
# 8    " Chi-squared Distance        ",
# 9    " Squared Chi-squared distance",
# 10   " Relative Euclidean Distance ",
# 11   " Relative Absolute Distance  ",
# 12   " Chord Distance              ",
# 13   " Squared Chord Distance      ",
# 14   " Geodesic Distance           ",
# 15   " Dice Coefficient            ",
# 16   " Jaccard Coefficient         ",
# 17   " Ochiai Coefficient          ",
# 18   " Edwards Chord               ",
# 19   " Overpeck Sq Chord Dist      ",


# TO DO - Fix problem calling paldist with integer data.

paldist <- function(y, dist.method="sq.chord") {
   METHODS <- c("euclidean", "sq.euclidean", "chord", "sq.chord", "chord.t", "sq.chord.t", "chi.squared", "sq.chi.squared", "bray")
   INDEX <- c(1, 2, 20, 19, 12, 13, 8, 9, 6)
   coef <- pmatch(dist.method, METHODS)
   if (is.na(coef)) {
      stop("Unknown distance measure")
   } else {
      coef <- INDEX[coef]
   }
   xx <- t(as.matrix(y))
   nr <- nrow(xx)
   nc <- ncol(xx)
   res <- as.double(matrix(NA, nrow=nc, ncol=nc))
   res <- .C("Dissim", as.vector(xx), as.vector(res), as.integer(nr), as.integer(nc), as.integer(coef), NAOK=TRUE, PACKAGE="rioja")[[2]]
   d <- as.dist(matrix(res, nrow = nc), diag=TRUE, upper=TRUE)
   attr(d, "Labels") <- colnames(xx)
   return(d)
}

paldist2 <- function(y1, y2, dist.method="sq.chord") {
   METHODS <- c("euclidean", "sq.euclidean", "chord", "sq.chord", "chord.t", "sq.chord.t", "chi.squared", "sq.chi.squared", "bray")
   INDEX <- c(1, 2, 20, 19, 12, 13, 8, 9, 6)
   coef <- pmatch(dist.method, METHODS)
   if (is.na(coef)) {
      stop("Unknown distance measure")
   } else {
      coef <- INDEX[coef]
   }
    xx1 <- as.matrix(t(y1))
    xx2 <- as.matrix(t(y2))
    nr1 <- nrow(xx1)
    nr2 <- nrow(xx2)
    nc1 <- ncol(xx1)
    nc2 <- ncol(xx2)
    if (nr1 != nr2)
        stop("number of columns different in y1 and y2")
    res <- as.double(matrix(NA, nrow = nc1, ncol = nc2))
    matrix(.C("Dissim2", as.vector(xx1), as.vector(xx2), as.vector(res), as.integer(nr1), as.integer(nc1), as.integer(nc2), as.integer(coef), NAOK=TRUE, PACKAGE="rioja")[[3]], nrow = nc1)
}
