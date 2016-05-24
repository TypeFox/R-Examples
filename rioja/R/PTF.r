# Utility functions for Quant Reconstructions II

rand.t.test <- function(object, ...) UseMethod("rand.t.test")

rand.t.test.default <- function(object, ...) { 
  stop(paste("No method to rand.t.test an object of class", class(object)[1]))
}

crossval <- function(object, ...) UseMethod("crossval")

crossval.default <- function(object, ...) {
  stop(paste("No method to cross-validate an object of class", class(object)[1]))
}

performance <- function(object, ...) UseMethod("performance")

performance.default <- function(object, ...) {
  stop(paste("No method for calculating performance for object of class", class(object)[1]))
}

.check.data <- function(y, x) {
  if (any(apply(y, 1, sum) < 1.0E-8))
     stop(paste("Species data have zero abundaces for the following rows:", which(apply(y, 1, sum) < 1.0E-8)))
  if (any(apply(y, 2, sum) < 1.0E-8))
     stop(paste("Species data have zero abundaces for the following columns:", which(apply(y, 2, sum) < 1.0E-8)))
}

.print.crossval <- function(object) {
  if (object$cv.summary$cv.method == "none" | object$cv.summary$cv.method == "loo")
    cat(paste("Cross val.         :", object$cv.summary$cv.method, "\n\n"))
  else if (object$cv.summary$cv.method == "lgo")
    cat(paste("Cross val.         :", object$cv.summary$cv.method, ": no. groups = ", object$cv.summary$ngroups , "\n\n"))
  else if (object$cv.summary$cv.method == "bootstrap")
    cat(paste("Cross val.         :", object$cv.summary$cv.method, ": no. boot cycles = ", object$cv.summary$nboot , "\n\n"))
  else if (object$cv.summary$cv.method == "h-block")
    cat(paste("Cross val.         :", object$cv.summary$cv.method, ": h-block cutoff = ", object$cv.summary$h.cutoff , "\n\n"))
}

.print.performance <- function(object, ...) {
   perf <- performance(object)
   results <- perf$object
   if (object$cv.summary$cv.method != "none") {
      results <- rbind(results, perf$crossval)
      n <- nrow(results)
      rownames(results)[(n/2+1):n] <- paste(rownames(results)[(n/2+1):n], "_XVal", sep="")
   }
   print.default(round(results, 4), print.gap = 2, ...)
}

.performance <- function(object, ...) {
  RMSE <- apply(residuals(object), 2, .rmse)
  R2 <- apply(object$fitted.values, 2, .r2, x=object$x)
  Avg.Bias <- apply(residuals(object), 2, mean, na.rm=TRUE)
  Max.Bias <- apply(residuals(object), 2, .max.bias, x=object$x)
  Skill <- 100 * (1.0 - apply((object$x - object$fitted)^2, 2, sum) / sum((object$x - mean(object$x))^2) )
  RMSE0 <- sqrt(mean((mean(object$x)-object$x)^2))
  res <- cbind(RMSE, R2, Avg.Bias, Max.Bias, Skill)
  result <- list(RMSE0=RMSE0, object=res)
  if (object$cv.summary$cv.method != "none") {
    if (object$cv.summary$cv.method == "bootstrap")
      RMSE <- object$cv.summary$RMSE.boot
    else
      RMSE <- apply(object$x - object$predicted, 2, .rmse)
    R2 <- apply(object$predicted, 2, .r2, x=object$x)
    Avg.Bias <- apply(object$x - object$predicted, 2, mean, na.rm=TRUE)
    Max.Bias <- apply(object$x - object$predicted, 2, .max.bias, x=object$x)
    Skill <- 100 * (1.0 - apply((object$x - object$predicted)^2, 2, sum) / sum((object$x - mean(object$x))^2))
    result.cv <- cbind(RMSE, R2, Avg.Bias, Max.Bias, Skill)
    result$crossval <- result.cv
  }
  result
}

.predict <- function(object, newdata=NULL, sse=FALSE, nboot=100, match.data=TRUE, verbose=TRUE, ...) {
  if (is.null(newdata)) {
     return (object$fitted.values)
  }
  if (is.null(newdata))
     return(object$fitted.values)
  if (match.data) {
    nms <- rownames(coef(object))
    mt <- match(colnames(newdata), nms)
    mt1 <- na.omit(mt)
    if (length(mt1) == 0)
       stop("No species in common between model and newdata")
    d <- matrix(0, nrow=nrow(newdata), ncol=length(nms))
    d[, mt1] <- as.matrix(newdata)[, !is.na(mt)]
    rownames(d) <- rownames(newdata)
    newdata <- d
  } else {
    if (ncol(object$y) != ncol(newdata)) 
       stop("Number of taxa does not match between datasets")
    if (any(colnames(object$y) != colnames(newdata)))
       stop("Taxon names do not match between datasets")
  }
#  nm.mod <- rownames(coef(object))
#  nm.new <- colnames(newdata)
#  mt <- match(nm.new, nm.mod)
#  mt1 <- na.omit(mt)
#  if (length(mt1) == 0 )
#    stop("Cannot predict: no taxa in common between object and newdata")
#  newdata2 <- matrix(0, ncol=nrow(coef(object)), nrow=nrow(newdata))
#  mt3 <- which(!is.na(mt))
#  newdata2[, mt1] <- as.matrix(newdata)[, mt3]
  predict.func <- paste("predict.internal", class(object)[1], sep=".")
  xHat.new <- do.call(predict.func, args=list(object=quote(object), y=quote(newdata), lean=FALSE, ...))
  rownames(xHat.new) <- rownames(newdata)
  colnames(xHat.new) <- colnames(object$fitted.values)
  if (sse) {
    if (verbose) {
      writeLines("Bootstrapping for SSE:")
      pb <- txtProgressBar(min = 0, max = 1, style = 3)
      on.exit(close(pb))
    }
    nsam <- nrow(object$y)
    nsam.new <- nrow(newdata)
    nest <- ncol(object$fitted.values)
    res2 <- array(dim=c(nsam, nest, nboot))
    res2.new <- array(dim=c(nsam.new, nest, nboot))
    call <- as.call(object$call.fit)
#    .set.rand.seed(100)
    for (i in 1:nboot) {
      o <- sample(nsam, replace=TRUE)
#      o <- apply(data.frame(rep(nsam, nsam)), 1, .get.rand) + 1
      out <- (1:nsam)[-unique(o)]
      y <- object$y[o, ]
      x <- object$x[o]
      y.test <- object$y[out, , drop=FALSE]
      mod <- eval(call)
      res2[out, , i] <- do.call(predict.func, args=list(object=quote(mod), y=quote(y.test), lean=TRUE, ...))
      res2.new[, , i] <- do.call(predict.func, args=list(object=quote(mod), y=quote(newdata), lean=TRUE, ...))
      if (verbose) {
        setTxtProgressBar(pb, i/nboot)
      }        
    }
    xHat <- object$fitted.values
    xHat.boot <- apply(res2, c(1,2), mean, na.rm=TRUE)
    xHat.new.boot <- apply(res2.new, c(1,2), mean, na.rm=TRUE)
    colnames(xHat.new.boot) <- colnames(xHat)
    rownames(xHat.new.boot) <- rownames(newdata)
    v1.boot <- apply(res2.new, c(1,2), sd, na.rm=TRUE)
    v2.boot <- apply(object$x-xHat.boot, 2, .rmse)
    colnames(v1.boot) <- colnames(xHat)
    rownames(v1.boot) <- rownames(newdata)
    SEP.boot <- sqrt(sweep(v1.boot^2, 2, v2.boot^2, "+"))
    colnames(SEP.boot) <- colnames(xHat)
    rownames(SEP.boot) <- rownames(newdata)
    results <- list(fit=xHat.new, fit.boot=xHat.new.boot, v1.boot=v1.boot, v2.boot=v2.boot, SEP.boot=SEP.boot)
  } else {
    results <- list(fit=xHat.new)
  }
  results
}

.crossval <- function(object, cv.method="loo", ngroups=10, nboot=100, verbose=TRUE, h.cutoff=0, h.dist, ...) 
{
  if (!("y" %in% names(object)))
    stop("Object does not contain species data, refit object using option lean=FALSE")
  METHODS <- c("loo", "lgo", "bootstrap", "h-block")
  cv.method <- pmatch(cv.method, METHODS)
  if (is.na(cv.method))
     stop("Unknown cross-validation method")
  nsam <- length(object$x)
  nres <- ncol(object$fitted.values)
#  func <- object$call.fit[1]
  call <- as.call(object$call.fit)
  predict.func <- paste("predict.internal", class(object)[1], sep=".")
  result <- matrix(nrow=nsam, ncol=nres)
  object$cv.summary$cv.method=METHODS[cv.method]
  if (verbose) {
    writeLines("Cross-validating:")
    pb <- txtProgressBar(min = 0, max = 1, style = 3)
    on.exit(close(pb))
  }
  if (cv.method == 1) {
    for (i in 1:nsam) {
      y <- object$y[-i, ]
      x <- object$x[-i]
      y.test <- object$y[i, , drop=FALSE]
      mod <- eval(call)
      xHat <- do.call(predict.func, args=list(object=quote(mod), y=quote(y.test), lean=TRUE, ...))
#      xHat <- do.call(predict, args=list(object=quote(mod), y=quote(y.test), lean=TRUE))
      result[i, ] <- xHat
      if (verbose) {
        setTxtProgressBar(pb, i/nsam)
      }        
    }
  } else if (cv.method == 2) {
    if (length(ngroups) > 1) {
      if (length(ngroups) != nsam)
        stop("length of leave-out groups does not equal number of samples")
      ngroups <- factor(ngroups)
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
      y <- object$y[-out, ]
      x <- object$x[-out]
      y.test <- object$y[out, , drop=FALSE]
      mod <- eval(call)
      xHat <- do.call(predict.func, args=list(object=quote(mod), y=quote(y.test), lean=TRUE))
      result[out, ] <- xHat
      if (verbose) {
        setTxtProgressBar(pb, i/ngroups)
      }
    }
    object$cv.summary$ngroups=ngroups
  } else if (cv.method == 3) {
    nest <- ncol(object$fitted.values)
    res2 <- array(dim=c(nsam, nest, nboot))
#    .set.rand.seed(100)
    for (i in 1:nboot) {
#      o <- apply(data.frame(rep(nsam, nsam)), 1, .get.rand) + 1
      o <- sample(nsam, replace=TRUE)
      out <- (1:nsam)[-unique(o)]
      y <- object$y[o, ]
      x <- object$x[o]
      y.test <- object$y[out, , drop=FALSE]
      mod <- eval(call)
      xHat <- do.call(predict.func, args=list(object=quote(mod), y=quote(y.test), lean=TRUE))
      res2[out, , i] <- xHat
      if (verbose) {
        setTxtProgressBar(pb, i/nboot)
      }
    }
    result <- apply(res2, c(1,2), mean, na.rm=TRUE)
    MS <- apply((object$x-res2)^2, c(1,2), mean, na.rm=TRUE)
    RMSE.boot <- sqrt(apply(MS, 2, mean, na.rm=TRUE))
    object$cv.summary$nboot=nboot
    object$cv.summary$RMSE.boot <- RMSE.boot
  } else if (cv.method == 4) {
    if (is.null(h.dist))
       stop("h-block cross-validation requested but h.dist is null") 
    h.dist <- as.matrix(h.dist)  
    if (nrow(h.dist) != ncol(h.dist))
       stop("h.dist doesn't look like a matrix of inter-site distances") 
    if (nrow(h.dist) != nsam) 
       stop(paste("Number of rows in h.dist (", nrow(h.dist), ") not equal to number of samples (", nsam, ")", sep="")) 
    nSamp <- vector("numeric", length=nsam)
    for (i in 1:nsam) {
       d <- h.dist[i, ]  
       sel <- d > h.cutoff
       y <- object$y[sel, , drop=FALSE]
       nSamp[i] <- nrow(y)
       if (nSamp[i] > 0) {
          keep <- abs(apply(y, 2, max)) > 0
          y <- y[, keep]
          x <- object$x[sel]
          y.test <- object$y[i, keep, drop=FALSE]
          mod <- eval(call)
          xHat <- do.call(predict.func, args=list(object=quote(mod), y=quote(y.test), lean=TRUE, ...))
          result[i, ] <- xHat
       } else {
          nmiss <- nmiss + 1
       }
       if (verbose) {
         setTxtProgressBar(pb, i/nsam)
       }
    }
    if (sum(nSamp < 1) > 0) {
       warning(paste(nmiss, "samples had no training samples with distance greater than ", h.cutoff, " and have not been predicted"))
    }
    object$n.h.block <- nSamp
    object$cv.summary$h.cutoff=h.cutoff
  } 
  colnames(result) <- colnames(object$fitted.values)
  rownames(result) <- rownames(object$fitted.values)
  object$predicted <- result
  object$residuals.cv <- object$x - result
  object
}
