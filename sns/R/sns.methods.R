######################################
#                                    #
#  Methods for sns objects           #
#                                    #
######################################

# convenience function for partitioning the state space
sns.make.part <- function(K, nsubset, method = "naive") {
  if (method != "naive") stop("invalid method") #TODO: consider implementing more sophisticated partitioning methods
  if (nsubset > K) stop("number of partitions cannot exceed state space dimensionality")
  
  # deterimining number of coordinates per subset (nvec)
  nvec <- rep(0, nsubset)
  nleft <- K
  c <- 1
  while (nleft > 0) {
    nvec[c] <- nvec[c] + 1
    nleft <- nleft - 1
    c <- c %% nsubset + 1
  }

  # assigning coordinates to subsets
  ret <- list()
  c <- 0
  for (n in 1:nsubset) {
    ret[[n]] <- as.integer(c + 1:nvec[n])
    c <- c + nvec[n]
  }
  if (sns.check.part(ret, K)) return (ret)
  else stop("unexpectedly invalid state space partitioning")
}

# function for checking that state space partitioning is valid
# (mutually-exclusive and collectively-exhaustive)
sns.check.part <- function(part, K) {
  return (identical(as.integer(sort(unlist(part))), 1:K))
}

# validating twice-differentiability and concavity of log-density
sns.check.logdensity <- function(x, fghEval
  #, numderiv = 0
  , numderiv.method = c("Richardson", "complex"), numderiv.args = list()
  , blocks = append(list(1:length(x)), as.list(1:length(x)))
  , dx = rep(1.0, length(x))
  , nevals = 100, negdef.tol = 1e-8
  , ...) {
  dx <- rep(dx, length(x)) # in case scalar dx is upplied by user
  fgh.ini <- fghEval(x, ...)
    
  ret <- list(check.ld.struct = NA, numderiv = NA, check.length.g = NA, check.dim.h = NA
              , x.mat = NA, fgh.ini = NA, fgh.num.ini = NA
              #, fgh.list = NA, fgh.num.list = NA
              , f.vec = NA, g.mat = NA, g.mat.num = NA, h.array = NA, h.array.num = NA
              , t.evals = NA, t.num.eval = NA
              , is.g.num.finite = NA, is.h.num.finite = NA
              , is.g.finite = NA, is.h.finite = NA
              , g.diff.max = NA, h.diff.max = NA
              , is.negdef.num = NA, is.negdef = NA
              )
  attr(ret, "class") <- "sns.check.logdensity"
  
  # determine if analytical derivatives are provided
  # numderiv = 0: g,h analytical
  # numderiv = 1: g analytical
  # numderiv = 2: no analytical
  is.fgh.list <- is.list(fgh.ini)
  if (is.fgh.list) {
    fgh.list.names <- names(fgh.ini)
    have.fgh <- c("f", "g", "h") %in% fgh.list.names
    if (all(have.fgh)) { # no numerical differentiation needed
      ret$check.ld.struct <- TRUE
      ret$numderiv <- 0
    } else if (all(have.fgh[1:2])) { # numerical gradient needed
      ret$check.ld.struct <- TRUE
      ret$numderiv <- 1
    } else {
      ret$check.ld.struct <- FALSE
      return (ret)
    }
  } else { # need numerical gradient and Hessian
    ret$check.ld.struct <- TRUE
    ret$numderiv <- 2
  }
  
  # fghEval with numerical differentiation
  numderiv.method <- match.arg(numderiv.method)
  if (ret$numderiv == 2) {
    fghEval.num <- function(x, ...) {
      f <- fghEval(x, ...)
      g <- grad(func = fghEval, x = x, ..., method = numderiv.method, method.args = numderiv.args)
      h <- hessian(func = fghEval, x = x, ..., method = numderiv.method, method.args = numderiv.args)
      return (list(f = f, g = g, h = h))
    }
  } else {
    fghEval.num <- function(x, ...) {
      f <- fghEval(x, ...)$f
      g <- grad(func = function(x, ...) fghEval(x, ...)$f, x = x, ..., method = numderiv.method, method.args = numderiv.args)
      h <- hessian(func = function(x, ...) fghEval(x, ...)$f, x = x, ..., method = numderiv.method, method.args = numderiv.args)
      return (list(f = f, g = g, h = h))
    }
  }
  
  fgh.num.ini <- fghEval.num(x, ...)
  
  #return (fgh.num.ini)
  
  # check dimensional consistency of gradient and Hessian
  K <- length(x)
  if (ret$numderiv == 0) {
    ret$check.length.g <- length(fgh.ini$g) == K
    ret$check.dim.h <- nrow(fgh.ini$h) == K && ncol(fgh.ini$h) == K
  } else if (ret$numderiv == 1) {
    ret$check.length.g <- length(fgh.ini$g) == K
    ret$check.dim.h <- nrow(fgh.num.ini$h) == K && ncol(fgh.num.ini$h) == K
  } else if (ret$numderiv == 2) {
    ret$check.length.g <- length(fgh.num.ini$g) == K
    ret$check.dim.h <- nrow(fgh.num.ini$h) == K && ncol(fgh.num.ini$h) == K
  }
  
  ret$x.mat <- sapply(1:K, function(k) runif(nevals, min = x[k] - 0.5 * dx[k], max = x[k] + 0.5 * dx[k]))
  ret$x.mat[1, ] <- x
  
  t <- proc.time()[3]
  fgh.list <- lapply(1:nevals, function(n) fghEval(ret$x.mat[n, ], ...))
  ret$t.evals <- proc.time()[3] - t
  
  t <- proc.time()[3]
  fgh.num.list <- lapply(1:nevals, function(n) fghEval.num(ret$x.mat[n, ], ...))
  ret$t.num.evals <- proc.time()[3] - t
  
  # twice-differentiability
  f.vec <- sapply(fgh.num.list, function(u) u$f)
  index.f.finite <- which(is.finite(f.vec))
  n.f.finite <- length(index.f.finite)
  ret$f.vec <- f.vec
  
  g.mat.num <- t(sapply(fgh.num.list, function(u) u$g))
  ret$g.mat.num <- g.mat.num
  ret$is.g.num.finite <- all(is.finite(g.mat.num[index.f.finite, ]))
  
  h.array.num <- array(t(sapply(fgh.num.list, function(u) u$h)), dim = c(nevals, K, K))
  ret$h.array.num <- h.array.num
  ret$is.h.num.finite <- all(is.finite(h.array.num[index.f.finite, , ]))
  
  # closeness of analytical and numerical results (we need better metrics of difference)
  if (ret$numderiv < 2) {
    g.mat <- t(sapply(fgh.list, function(u) u$g))
    ret$g.mat <- g.mat
    ret$is.g.finite <- all(is.finite(g.mat[index.f.finite, ]))
    ret$g.diff.max <- max(sapply(1:n.f.finite, function(i)
      sqrt(sum((g.mat[index.f.finite[i], ] - g.mat.num[index.f.finite[i], ])^2))/sqrt(sum(g.mat[index.f.finite[i], ]^2))))
  }
  if (ret$numderiv == 0) {
    h.array <- array(t(sapply(fgh.list, function(u) u$h)), dim = c(nevals, K, K))
    ret$h.array <- h.array
    ret$is.h.finite <- all(is.finite(h.array[index.f.finite, , ]))
    ret$h.diff.max <- max(sapply(1:n.f.finite, function(i)
      norm(h.array[index.f.finite[i], , ] - h.array.num[index.f.finite[i], , ], type = "F")/norm(h.array[index.f.finite[i], , ], type = "F")))
  }
  
  # concavity (negative definiteness)
  sns.is.positive.definite <- function(A, tol = 1e-8) all(eigen(A)$values > tol)
  
  ret$is.negdef.num <- sapply(1:length(blocks), function(i) all(sapply(1:n.f.finite, function(u) {
    h.tmp <- h.array.num[index.f.finite[u], , ]
    sns.is.positive.definite(-h.tmp[blocks[[i]], blocks[[i]], drop = F], tol = negdef.tol)
  })))
  
  if (ret$numderiv == 0) {
    ret$is.negdef <- sapply(1:length(blocks), function(i) all(sapply(1:n.f.finite, function(u) {
      h.tmp <- h.array[index.f.finite[u], , ]
      sns.is.positive.definite(-h.tmp[blocks[[i]], blocks[[i]], drop = F], tol = negdef.tol)
    })))
  }
  
  return (ret)
}

print.sns.check.logdensity <- function(x, ...) {
  sns.TF.to.YesNo <- function(x) ifelse(x, "Yes", "No")
  if (!x$check.ld.struct) {
    cat("log-density output list has invalid names\n")
    return (invisible(NULL))
  }
  numderiv <- x$numderiv
  #cat("log-density is valid\n")
  cat("number of finite function evals:", length(which(is.finite(x$f.vec))), "(out of ", length(x$f.vec), ")\n")
  cat("recommended numderiv value:", numderiv, "\n")
  cat("is length of gradient vector correct?", sns.TF.to.YesNo(x$check.length.g), "\n")
  cat("are dims of Hessian matrix correct?", sns.TF.to.YesNo(x$check.dim.h), "\n")
  cat("is numerical gradient finite?", sns.TF.to.YesNo(x$is.g.num.finite), "\n")
  cat("is numerical Hessian finite?", sns.TF.to.YesNo(x$is.h.num.finite), "\n")
  if (numderiv < 2) {
    cat("is analytical gradient finite?", sns.TF.to.YesNo(x$is.g.finite), "\n")
    cat("maximum relative diff in gradient:", x$g.diff.max, "\n")
  }
  if (numderiv == 0) {
    cat("is analytical Hessian finite?", sns.TF.to.YesNo(x$is.h.finite), "\n")
    cat("maximum relative diff in Hessian:", x$h.diff.max, "\n")
  }
  cat("is numeric Hessian (block) negative-definite?", sns.TF.to.YesNo(x$is.negdef.num), "\n")
  if (numderiv == 0) {
    cat("is analytical Hessian (block) negative-definite?", sns.TF.to.YesNo(x$is.negdef), "\n")
  }
  
  return (invisible(NULL))
}

# predict methods
predict.sns <- function(object, fpred
  , nburnin = max(nrow(object)/2, attr(object, "nnr"))
  , end = nrow(object), thin = 1, ...) {

  niter <- nrow(object)
  nnr <- attr(object, "nnr")
  nmcmc <- niter - nnr
  if (nburnin < nnr) warning("it is strongly suggested that burnin period includes NR iterations (which are not valid MCMC iterations)")
  myseq <- seq(from = nburnin + 1, to = end, by = thin)

  pred <- apply(object[myseq, ], 1, fpred, ...)
  class(pred) <- "predict.sns"
  return (pred)
}
summary.predict.sns <- function(object, quantiles = c(0.025, 0.5, 0.975)
  , ess.method = c("coda", "ise"), ...) {
  smp.mean <- rowMeans(object)
  smp.sd <- apply(object, 1, sd)
  smp.ess <- ess(t(object), method = ess.method[1])
  smp.quantiles <- t(apply(object, 1, quantile, probs = quantiles))
  ret <- list(mean = smp.mean, sd = smp.sd, ess = smp.ess, quantiles = smp.quantiles, nseq = ncol(object))
  class(ret) <- "summary.predict.sns"
  return (ret)
}
print.summary.predict.sns <- function(x, ...) {
  cat("prediction sample statistics:\n")
  cat("\t(nominal sample size: ", x$nseq, ")\n", sep="")
  stats <- cbind(x$mean, x$sd, x$ess, x$quantiles)
  colnames(stats)[1:3] <- c("mean", "sd", "ess")
  rownames(stats) <- c(1:length(x$mean))
  printCoefmat(stats[1:min(length(x$mean), 6), ])
  if (length(x$mean) > 6) cat("...\n")
}

# print method 
print.sns <- function(x, ...) {
  cat("Stochastic Newton Sampler (SNS)\n")
  cat("state space dimensionality: ", ncol(x), "\n")
  if (!is.null(attr(x, "part"))) cat("state space partitioning: ", attr(x, "part"), " subsets\n")
  cat("total iterations: ", nrow(x), "\n")
  cat("\t(initial) NR iterations:", attr(x, "nnr"), "\n")
  cat("\t(final) MCMC iterations:", nrow(x) - attr(x, "nnr"), "\n")
}

# summary methods
# primary output:
# 1) acceptance rate
# 2) mean relative deviation (if available)
# 3) sample statistics (mean, sd, quantiles, ess, pval) (if available)
summary.sns <- function(object, quantiles = c(0.025, 0.5, 0.975)
  , pval.ref = 0.0, nburnin = max(nrow(object)/2, attr(object, "nnr"))
  , end = nrow(object), thin = 1, ess.method = c("coda", "ise"), ...) {
  K <- ncol(object)
  nnr <- attr(object, "nnr")
  if (nburnin < nnr) warning("it is strongly suggested that burnin period includes NR iterations (which are not valid MCMC iterations)")
  
  # number of subsets in state space partitioning
  npart <- max(1, length(attr(object, "part")))
    
  # average relative deviation of function value from quadratic approximation (post-burnin)
  if (!is.null(attr(object, "reldev"))) reldev.mean <- mean(attr(object, "reldev"), na.rm = TRUE)
  else reldev.mean <- NA
  
  nsmp <- end - nburnin
  if (nsmp > 0) {
    # average acceptance rate for MH transition proposals
    accept.rate <- sum(attr(object, "accept")[nburnin + 1:nsmp, ]) / length(attr(object, "accept")[nburnin + 1:nsmp, ])
    
    myseq <- seq(from = nburnin + 1, to = end, by = thin)
    nseq <- length(myseq)
    
    smp.mean <- colMeans(object[myseq, ])
    smp.sd <- apply(object[myseq, ], 2, sd)
    smp.ess <- ess(object[myseq, ], method = ess.method[1])
    smp.quantiles <- t(apply(object[myseq, ], 2, quantile, probs = quantiles))
    smp.pval <- apply(object[myseq, ], 2, sns.calc.pval, ref = pval.ref, na.rm = FALSE)
    
  } else {
    accept.rate <- NA
    nseq <- 0
    
    smp.mean <- NA
    smp.sd <- NA
    smp.ess <- NA
    smp.quantiles <- NA
    smp.pval <- NA
  }
  ret <- list(K = K, nnr = nnr, nburnin = nburnin, end = end, thin = thin
    , niter = nrow(object), nsmp = nsmp, nseq = nseq, npart = npart
    , accept.rate = accept.rate, reldev.mean = reldev.mean
    , pval.ref = pval.ref, ess.method = ess.method
    , smp = list(mean = smp.mean, sd = smp.sd, ess = smp.ess, quantiles = smp.quantiles, pval = smp.pval))
  class(ret) <- "summary.sns"
  return (ret)
}

print.summary.sns <- function(x, ...) {
  cat("Stochastic Newton Sampler (SNS)\n")
  cat("state space dimensionality: ", x$K, "\n")
  if (x$npart > 1) cat("state space partitioning: ", x$npart, " subsets\n")
  cat("total iterations: ", x$niter, "\n")
  cat("\tNR iterations: ", x$nnr, "\n")
  cat("\tburn-in iterations: ", x$nburnin, "\n")
  cat("\tend iteration: ", x$end, "\n")
  cat("\tthinning interval: ", x$thin, "\n")
  cat("\tsampling iterations (before thinning): ", x$nsmp, "\n")
  #cat("\tsampling iterations (after thinning): ", x$nseq, "\n")
  cat("acceptance rate: ", x$accept.rate, "\n")
  if (!is.na(x$reldev.mean)) cat("\tmean relative deviation from quadratic approx:", format(100*x$reldev.mean, digits=3), "% (post-burnin)\n")
  if (x$nsmp > 0) {
    cat("sample statistics:\n")
    cat("\t(nominal sample size: ", x$nseq, ")\n", sep="")
    stats <- cbind(x$smp$mean, x$smp$sd, x$smp$ess, x$smp$quantiles, x$smp$pval)
    colnames(stats)[c(1:3, 4 + ncol(x$smp$quantiles))] <- c("mean", "sd", "ess", "p-val")
    rownames(stats) <- c(1:x$K)
    printCoefmat(stats[1:min(x$K, 6), ], P.values = TRUE, has.Pvalue = TRUE)
    if (x$K > 6) cat("...\n")
    cat("summary of ess:\n")
    print(summary(x$smp$ess))
  }
}

# plot method
plot.sns <- function(x, nburnin = max(nrow(x)/2, attr(x, "nnr"))
  , select = if (length(x) <= 10) 1:5 else 1, ...) {
  init <- attr(x, "init")
  lp.init <- attr(x, "lp.init")
  lp <- attr(x, "lp")
  
  # in all cases, vertical line delineates transition from nr to mcmc mode
  K <- ncol(x)
  niter <- nrow(x)
  nnr <- attr(x, "nnr")
  if (nburnin < nnr) warning("it is strongly suggested that burnin period includes NR iterations (which are not valid MCMC iterations)")
  
  # log-probability trace plot
  if (1 %in% select) {
    plot(0:niter, c(lp.init, lp), type = "l"
      , xlab = "iter", ylab = "log-probability", main = "Log-Probability Trace Plot")
    if (nnr > 0 && nnr < niter) abline(v = nnr + 0.5, lty = 2, col = "red")
  }
  
  # state vector trace plots
  if (2 %in% select) {
    for (k in 1:K) {
      plot(0:niter, c(init[k], x[, k]), type = "l"
        , xlab = "iter", ylab = paste("x[", k, "]", sep = ""), main = "State Variable Trace Plot")
      if (nnr > 0 && nnr < niter) abline(v = nnr + 0.5, lty = 2, col = "red")
    }
  }
  
  if (nburnin < niter) {

    if (3 %in% select) {
      # effective sample size (horizontal line is maximum possible effective sample size)
      my.ess <- ess(x[(nburnin + 1):niter, ])
      plot(1:K, my.ess, xlab = "k", ylab = "effective sample size", ylim = c(0, niter - nburnin), main = "Effective Sample Size by Coordinate")
      abline(h = niter - nburnin, lty = 2, col = "red")
    }
    
    if (4 %in% select) {
      # state vector (univariate) histograms
      K <- ncol(x)
      for (k in 1:K) {
        hist(x[(nburnin + 1):niter, k], xlab = paste("x[", k, "]", sep = ""), main = "State Variable Histogram (post-burnin)")
        abline(v = mean(x[(nburnin + 1):niter, k]), lty = 2, col = "red")
      }
    }
  
    if (5 %in% select) {
      # state vector (univariate) autocorrelation plots
      K <- ncol(x)
      for (k in 1:K) {
        acf(x[(nburnin + 1):niter, k], xlab = paste("x[", k, "]", sep = ""), main = "State Variable Autocorrelation Plot (post-burnin)")
      }
    }

  }
}

sns.calc.pval <- function(x, ref=0.0, na.rm = FALSE) { # add flag for one-sided vs. two-sided
  if (na.rm) x <- x[!is.na(x)]
  bigger <- median(x)>ref
  if (sd(x)<.Machine$double.eps) {
    ret <- NA
  } else {
    ret <- max(1/length(x), 2*length(which(if (bigger) x<ref else x>ref))/length(x)) # TODO: justify minimum value
  }
  attr(ret, "bigger") <- bigger
  return (ret)
}

# convenience function for numerical augmentation of a log-density
sns.fghEval.numaug <- function(fghEval, numderiv = 0
  , numderiv.method = c("Richardson", "simple"), numderiv.args = list()) {
  numderiv <- as.integer(numderiv)
  if (numderiv > 0) {
    numderiv.method <- match.arg(numderiv.method)
    if (numderiv == 1) { # we need numeric hessian
      fghEval.int <- function(x, ...) {
        fg <- fghEval(x, ...)
        h <- hessian(func = function(x, ...) fghEval(x, ...)$f, x = x, ..., method = numderiv.method, method.args = numderiv.args)
        return (list(f = fg$f, g = fg$g, h = h))
      }
    } else { # we need numeric gradient and hessian
      fghEval.int <- function(x, ...) {
        f <- fghEval(x, ...)
        g <- grad(func = fghEval, x = x, ..., method = numderiv.method, method.args = numderiv.args)
        h <- hessian(func = fghEval, x = x, ..., method = numderiv.method, method.args = numderiv.args)
        return (list(f = f, g = g, h = h))
      }
    }
  } else {
    fghEval.int <- fghEval
  }
}


