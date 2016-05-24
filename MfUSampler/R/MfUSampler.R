# internal function: wrapper around a multivariate function to convert it to univariate, to be used with univariate slice sampler
MfU.fEval <- function(xk, k, x, f, ...) {
  x[k] <- xk
  return (f(x, ...))
}

# internal functions: wrappers around a multvariate function to extract function and gradient, to be used with adaptive rejection sampler
# we need to implement a vectorized version, since ARS code expects it
MfU.fgEval.f <- function(uk, k, u, func, ...) {
  ret <- sapply(uk, function(ukk) {
    u[k] <- ukk
    return (func(u, ..., grad=FALSE))
  })
  return (ret)
}
MfU.fgEval.g <- function(uk, k, u, func, ...) {
  ret <- sapply(uk, function(ukk) {
    u[k] <- ukk
    return (func(u, ..., grad=TRUE)[k])
  })
  return (ret)
}

# public function: multivariate sampler, utilizing univariate samplers (slice/ars/arms) through a Gibbs wrapper
MfU.Sample <- function(x, f, uni.sampler = c("slice", "ars", "arms", "unimet")
  , ..., control = MfU.Control(length(x))) {
  if (uni.sampler == "slice") {
    for (k in 1:length(x)) {
      x[k] <- MfU.UniSlice(x[k], MfU.fEval, k, x, f, ..., w=control$slice$w[k], m=control$slice$m[k]
                           , lower=control$slice$lower[k], upper=control$slice$upper[k])
    }
    return (x)
  } else if (uni.sampler == "ars") {
    for (k in 1:length(x)) {
      x[k] <- ars(n=1, MfU.fgEval.f, MfU.fgEval.g, x=control$ars$x[[k]], ns=control$ars$ns[k], m=control$ars$m[k]
                  , emax=control$ars$emax[k], lb=control$ars$lb[k], ub=control$ars$ub[k], xlb=control$ars$xlb[k]
                  , xub=control$ars$xub[k], k, x, f, ...)
    }
    return (x)
  } else if (uni.sampler == "arms") {
    for (k in 1:length(x)) {
      x[k] <- arms(x[k], MfU.fEval, indFunc = control$arms$indFunc[[k]]
                   , n.sample = 1, k, x, f, ...)
    }
    return (x)
  } else if (uni.sampler == "unimet") {
    for (k in 1:length(x)) {
      x[k] <- MfU.UniMet(x[k], MfU.fEval, k, x, f, ..., sigma = control$unimet$sigma[k])
    }
    return (x)
  } else {
    stop("invalid sampler")
  }
}

# public function: convenience wrapper around MfU.Sample to obtain multiple samples
MfU.Sample.Run <- function(x, f, uni.sampler = c("slice", "ars", "arms", "unimet")
  , ..., control = MfU.Control(length(x)), nsmp = 10) {
  uni.sampler <- match.arg(uni.sampler)
  
  t <- proc.time()[3]
  x.smp <- array(NA, dim = c(nsmp, length(x)))
  for (n in 1:nsmp) {
    x <- MfU.Sample(x = x, f = f, uni.sampler = uni.sampler, control = control, ...)
    x.smp[n, ] <- x
  }
  t <- proc.time()[3] - t
  
  class(x.smp) <- c("MfU", class(x.smp))
  attr(x.smp, "t") <- t
  return (x.smp)
}

# generic methods
# summary
summary.MfU <- function(object, start = round(nrow(object)/2) + 1, end = nrow(object), thin = 1
  , quantiles = c(0.025, 0.5, 0.975)
  , ...) {
  coda.object <- mcmc(data = as.matrix(object), start = start, end = end, thin = thin)
  ret <- summary(coda.object, quantiles = quantiles)
  
  # adding a few items to coda output
  myseq <- seq(from = start, to = end, by = thin)
  # sample covariance matrix
  ret$covar <- cov(object[myseq, ])
  # effective sample size
  ess <- effectiveSize(object[myseq, ])
  t <- attr(object, "t")
  # fraction of total sampling time assigned to selected samples
  t.myseq <- t * length(myseq) / nrow(object)
  # number of independent samples per sec
  iss <- ess / t.myseq
  
  ret$t <- t
  ret$t.myseq <- t.myseq
  ret$ess <- ess
  ret$iss <- iss
  
  ret$ntot <- nrow(object)
  ret$nseq <- length(myseq)
  
  class(ret) <- "summary.MfU"
  
  return (ret)
}
print.summary.MfU <- function(x, ...) {
  x.coda <- x
  class(x.coda) <- "summary.mcmc"
  print(x.coda)
  cat("time for all samples (", x$ntot, "):", x$t, "sec\n")
  cat("time assigned to selected samples (", x$nseq, "):", x$t.myseq, "sec\n")
  cat("Effective sample size / independent samples per sec:\n")
  df <- cbind(x$ess, x$iss)
  colnames(df) <- c("ess", "iss")
  print(df)
}

# plot
plot.MfU <- function(x, start = round(nrow(x)/2) + 1
  , end = nrow(x), thin = 1, ...) {
  x <- mcmc(data = as.array(x), start = start, end = end, thin = thin)
  plot(x, ...)
}

# predict
predict.MfU <- function(object, fpred, ...) {
  pred <- t(apply(object, 1, fpred, ...))
  class(pred) <- c("predict.MfU", class(pred))
  return (pred)
}
summary.predict.MfU <- function(object
  , start = round(nrow(object)/2) + 1
  , end = nrow(object), thin = 1
  , quantiles = c(0.025, 0.5, 0.975), ...) {
  
  myseq <- seq(from = start, to = end, by = thin)
  object <- object[myseq, ]
  
  smp.mean <- colMeans(object)
  smp.sd <- apply(object, 2, sd)
  smp.ess <- effectiveSize(object)
  smp.quantiles <- t(apply(object, 2, quantile, probs = quantiles))
  ret <- list(mean = smp.mean, sd = smp.sd, ess = smp.ess, quantiles = smp.quantiles, nseq = nrow(object))
  class(ret) <- "summary.predict.MfU"
  return (ret)
}
print.summary.predict.MfU <- function(x, n = 6L, ...) {
  cat("prediction sample statistics:\n")
  cat("\t(nominal sample size: ", x$nseq, ")\n", sep="")
  stats <- cbind(x$mean, x$sd, x$ess, x$quantiles)
  colnames(stats)[1:3] <- c("mean", "sd", "ess")
  rownames(stats) <- c(1:length(x$mean))
  printCoefmat(stats[1:min(length(x$mean), n), ])
  if (length(x$mean) > n) cat("...\n")
}

# public function: setting tuning parameters of univariate samplers
MfU.Control <- function(n, slice.w=1, slice.m=Inf, slice.lower=-Inf, slice.upper=+Inf
  , ars.x=c(-4,1,4), ars.ns=100, ars.m=3, ars.emax=64, ars.lb=FALSE, ars.xlb=0, ars.ub=FALSE, ars.xub=0
  , arms.indFunc = function(x) TRUE, unimet.sigma = 1.0) {
  
  if (missing(n)) stop("dimensionality of state space (n) must be provided")
  
  if (is.list(arms.indFunc)) {
    arms.indFunc.aug <- list()
    for (i in 1:length(arms.indFunc)) {
      arms.indFunc.aug[[i]] <- function(xk, k, x, f, ...) {
        x[k] <- xk
        return (arms.indFunc[[i]](x, ...))
      }
    }
  } else {
    arms.indFunc.aug <- function(xk, k, x, f, ...) {
      x[k] <- xk
      return (arms.indFunc(x, ...))
    }
  }
  
  expand.to.vector <- function(x, n) {
    if (length(x)==1) {
      return (rep(x,n))
    } else if (length(x)==n) {
      return (x)
    } else {
      stop("invalid x")
    }
  }
  expand.to.list <- function(x, n) {
    if (!is.list(x)) {
      ret <- list()
      for (i in 1:n) ret[[i]] <- x
      return (ret)
    } else if (length(x)==n) {
      return (x)
    } else {
      stop("invalid x")
    }
  }
  list(slice=list(w=expand.to.vector(slice.w, n), m=expand.to.vector(slice.m, n)
                  , lower=expand.to.vector(slice.lower, n), upper=expand.to.vector(slice.upper, n))
       , ars=list(x=expand.to.list(ars.x, n), ns=expand.to.vector(ars.ns, n), m=expand.to.vector(ars.m, n)
                  , emax=expand.to.vector(ars.emax, n), lb=expand.to.vector(ars.lb, n), ub=expand.to.vector(ars.ub, n)
                  , xlb=expand.to.vector(ars.xlb, n), xub=expand.to.vector(ars.xub, n))
       , arms = list(indFunc = expand.to.list(arms.indFunc.aug, n))
       , unimet = list(sigma = expand.to.vector(unimet.sigma, n)))
}
