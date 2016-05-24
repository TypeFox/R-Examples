LWR <- function(y, x, FUN=WA, dist.method="sq.chord", k=30, lean=TRUE, fit.model=TRUE, check.data=TRUE, verbose=TRUE, ...)
{
  call.fit <- as.call(list(quote(LWR), y=quote(y), x=quote(x), FUN=substitute(FUN), dist.method=dist.method, k=k, lean=lean))
  funname <- substitute(FUN)
  FUN <- match.fun(FUN) 
  if (check.data) {
    if (any(apply(y, 1, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following rows:", paste(which(apply(y, 1, sum) < 1.0E-8), collapse=",")))
    if (any(apply(y, 2, sum) < 1.0E-8))
       stop(paste("Species data have zero abundances for the following columns:", paste(which(apply(y, 2, sum) < 1.0E-8), collapse=",")))
  }
  if (k >= nrow(y))
    stop("k is greater than number of training samples")
  x1 <- as.numeric(x)
  n <- 2
  diss <- as.matrix(paldist(y, dist.method=dist.method))
  ind <- apply(diss, 2, order)
  ind <- ind[n:(n+k-1), ]
  dist.n <- t(apply(diss, 2, sort)[n:(n+k-1), ])
  rownames(dist.n) <- rownames(y)
  colnames(dist.n) <- paste("N", sprintf("%02d", 1:k), sep="")
  nsam <- nrow(y)
  ncol.res <- NA
  cut1 <- 0
  cut2 <- 0
  if (verbose) {
    writeLines("Fitting models:")
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    on.exit(close(pb))
  }
  
  if (funname=="MLRC") {
     cut1 <- 2
     cut2 <- 0.01
  }
  if (fit.model) {
    for (i in 1:nsam) {
      if (verbose) {
        setTxtProgressBar(pb, i/nsam)
      }
      ytrain <- y[ind[, i], ]
      N <- apply(ytrain>0, 2, sum)
      Mx <- apply(ytrain, 2, sum)
      xtrain <- x[ind[, i]]
      if (funname=="MLRC") {
         ytrain <- ytrain[, N>cut1 & Mx > cut2]
         mod <- FUN(ytrain, xtrain, check.data=FALSE, ...)
         pred <- predict(mod, y[i, N>cut1 & Mx>cut2, drop=FALSE], verbose=FALSE)$fit
      } else {
         mod <- FUN(ytrain, xtrain, check.data=FALSE, ...)
         pred <- predict(mod, y[i, , drop=FALSE], verbose=FALSE)$fit
      }
      if (i == 1) {
         ncol.res <- ncol(pred)
         xHat <- matrix(NA, nrow=nrow(y), ncol=ncol.res)
         colnames(xHat) <- colnames(pred)
      }
      xHat[i, ] <- pred
    }
  }
  rownames(xHat) <- rownames(y)
  nms <- t(matrix(rownames(y)[ind], nrow=k))
  rownames(nms) <- rownames(y)
  rownames(diss) <- rownames(y)
  colnames(diss) <- rownames(y)
  call.print <- match.call()
  result <- list(call=call.fit, call.print=call.print, fitted.values=xHat, dist.n=dist.n, ind=ind, match.name=nms, x=x1, dist.method=dist.method, k=k, y=y, FUN=FUN)
  if (!lean)
     result <- c(result, list(dist=diss))
  result$cv.summary <- list(cv.method="none")
	class(result) <- "LWR"
	result
}

predict.LWR <- function(object, newdata=NULL, k=object$k, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, lean=TRUE, ...)
{
  if (k < 1 | k > object$k)
    stop("k out of range")
  FUN <- match.fun(object$call$FUN)
  if (is.null(newdata))
     return(object$fitted.values)
  if (k >= nrow(object$y))
    stop("k is less than number of training samples")
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
  diss <- paldist2(y1, y2, dist.method=object$dist.method)
  ind <- apply(diss, 2, order)
  ind <- ind[n:(n+k-1), , drop=FALSE]
  dist.n <- t(apply(diss, 2, sort)[n:(n+k-1), , drop=FALSE])
  rownames(dist.n) <- rownames(y2)
  colnames(dist.n) <- paste("N", sprintf("%02d", 1:k), sep="")
  ncol.res <- NA
  for (i in 1:nrow(y2)) {
    ytrain <- y1[ind[, i], , drop=FALSE]
    xtrain <- x1[ind[, i]]
    mx <- apply(ytrain, 2, sum)
    ytrain <- ytrain[, mx>0]
    ytest <- y2[i, mx>0, drop=FALSE]
    mod <- FUN(ytrain, xtrain, check.data=FALSE, ...)
    pred <- predict(mod, ytest, verbose=FALSE)$fit
    if (i == 1) {
       ncol.res <- ncol(pred)
       xHat <- matrix(NA, nrow=nrow(y2), ncol=ncol.res)
       colnames(xHat) <- colnames(pred)
    }
    xHat[i, ] <- pred
  }
  rownames(xHat) <- rownames(y2)
  nms <- t(matrix(rownames(y1)[ind[n:(n+k-1), ]], nrow=k))
  rownames(nms) <- rownames(y2)
  colnames(nms) <- paste("N", sprintf("%02d", 1:k), sep="")
  rownames(diss) <- rownames(y1)
  colnames(diss) <- rownames(y2)
  result <- list(fit=xHat, dist.n=dist.n, match.name=nms)
  if (!lean)
     result <- c(result, list(dist=diss))
  if (sse) {
   
   stop("Sample Specific Errors not yet implemented for LWR")
     
   feedback <- ifelse(is.logical(verbose), 50, as.integer(verbose))
   if (is.null(object$dist))
      stop("No distances: refit original model using \"lean=FALSE\"")
    nsam <- nrow(object$y)
    nsam.new <- nrow(newdata)
    nest <- ncol(xHat)
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
      dist.n <- t(apply(diss.m, 2, sort)[1:k, ])
      x.n <- t(matrix(x[ind[1:k, ]], nrow=k))
      res2[out, 1, i] <- apply(x.n, 1, mean, na.rm=TRUE)
      res2[out, 2, i] <- rowSums((x.n / dist.n), na.rm=TRUE) / rowSums(1/dist.n)
      
      diss.f <- diss[o, ]
      ind <- apply(diss.f, 2, order)
      dist.n <- apply(diss.f, 2, sort)[1:k, ]
      x.n <- matrix(x[ind[1:k, ]], nrow=k)
      
# this is the code for MAT - it needs updating      

#      res2.new[, 1, i] <- apply(x.n, 2, mean, na.rm=TRUE)
#      res2.new[, 2, i] <- colSums((x.n / dist.n), na.rm=TRUE) / colSums(1/dist.n)
      
      if (verbose) {
          if (i %% feedback == 0) {
            cat (paste("Bootstrap sample", i, "\n"))
            flush.console()
          }
      }
    }      
    xHat.boot <- apply(res2, c(1,2), mean, na.rm=TRUE)
    xHat.new.boot <- apply(res2.new, c(1,2), mean, na.rm=TRUE)
#    colnames(xHat.new.boot) <- colnames(result$fit)
#    rownames(xHat.new.boot) <- rownames(newdata)
    v1.boot <- apply(res2.new, c(1,2), sd, na.rm=TRUE)
    v2.boot <- apply(object$x-xHat.boot, 2, .rmse)
#    colnames(v1.boot) <- colnames(result$fit)
#    rownames(v1.boot) <- rownames(newdata)
    SEP.boot <- sqrt(sweep(v1.boot^2, 2, v2.boot^2, "+"))
#    colnames(SEP.boot) <- colnames(result$fit)
#    rownames(SEP.boot) <- rownames(newdata)
    result <- c(result, list(fit.boot=xHat.new.boot, v1.boot=v1.boot, v2.boot=v2.boot, SEP.boot=SEP.boot))
  }   
  result
}


crossval.LWR <- function(object, k=object$k, cv.method="lgo", verbose=TRUE, ngroups=10, nboot=100, h.cutoff=0, h.dist=NULL, ...)
{
  if (k < 1 | k > object$k)
    stop("k out of range")
  METHODS <- c("lgo", "bootstrap")
  cv.method <- pmatch(cv.method, METHODS)
  FUN <- match.fun(object$call$FUN)
  if (is.na(cv.method))
     stop("Unknown cross-validation method")
  nsam <- length(object$x)
  nres <- ncol(object$fitted.values)
  result <- matrix(nrow=nsam, ncol=nres)
  object$cv.summary$cv.method=METHODS[cv.method]
  k <- object$k
  if(verbose) {
    writeLines("Cross-validating:")
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    on.exit(close(pb))
  }
  if (cv.method == 1) {
    if (length(ngroups) > 1) {
       if (length(ngroups) != nsam)
          stop("Length of leave-out groups does not equal number of samples")
       grps <- as.integer(ngroups)
       ngroups <- length(unique(ngroups))
       o <- 1:nsam
    }
    else {
      o <- sample(nsam)
      grps <- rep(1:ngroups, length.out=nsam) 
    }
    for (i in 1:ngroups) {
      out <- o[grps==i]
      x.test <- object$x[out, drop=FALSE]
      x.train <- object$x[-out, drop=FALSE]
      y.test <- object$y[out, , drop=FALSE]
      y.train <- object$y[-out, , drop=FALSE]
      sel <- apply(y.train, 2, sum) > 0
      y.train <- y.train[, sel]
#      y.test <- y.test[, sel]
      mod <- LWR(y.train, x.train, FUN, dist.method=object$dist.method, k=k, lean=TRUE, fit.model=TRUE, check.data=TRUE, verbose=FALSE, ...)
      xHat <- predict.LWR(mod, y.test, verbose=FALSE)
      result[out, ] <- xHat$fit
      if (verbose) {
        setTxtProgressBar(pb, i/ngroups)
      }
    }
    object$cv.summary$ngroups=ngroups
  } else if (cv.method == 2) {
    res2 <- array(dim=c(nsam, nres, nboot))
    for (i in 1:nboot) {
      o <- sample(nsam, replace=TRUE)
      out <- (1:nsam)[-unique(o)]
      x.test <- object$x[out, drop=FALSE]
      x.train <- object$x[-out, drop=FALSE]
      y.test <- object$y[out, , drop=FALSE]
      y.train <- object$y[-out, , drop=FALSE]
      sel <- apply(y.train, 2, sum) > 0
      y.train <- y.train[, sel]
      mod <- LWR(y.train, x.train, FUN, dist.method=object$dist.method, k=k, lean=TRUE, fit.model=TRUE, check.data=TRUE, verbose=FALSE, ...)
      xHat <- predict.LWR(mod, y.test, verbose=FALSE)
      res2[out, , ] <- xHat$fit
      if (verbose) {
        setTxtProgressBar(pb, i/nboot)
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
  object$cv.summary$cv.method = cv.method
  object
}


print.LWR <- function(x, ...) {
  cat("\n")
  cat("Method : Locally Weighted Regression\n")
  cat("Call   : ")
  print(x$call)
##  cat(paste(deparse(x$call), "\n\n"))
  cat(paste("Distance           :", x$dist.method, "\n"))
  cat(paste("No. samples        :", length(x$x), "\n"))
  cat(paste("No. species        :", ncol(x$y), "\n"))
  cat(paste("No. local          :", x$k, "\n"))
  
  .print.crossval(x)
  cat("\nPerformance:\n")
  .print.performance(x)
  cat("\n")
}

performance.LWR <- function(object, ...) {
  .performance(object, ...)
}

summary.LWR <- function(object, full=FALSE, ...)
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

fitted.LWR <- function(object, ...) {
  object$fitted.values
}

residuals.LWR <- function(object, cv=FALSE, ...) {
  if (cv == FALSE)
     return (object$x - object$fitted.values)
  else {
     if (object$cv.summary$cv.method == "none")
        stop("Object does not contain cross validation results")
     return (object$residuals.cv)
  }
}

