confint.lnre <- function (object, parm, level=0.95, method=c("normal", "mad", "empirical"), plot=FALSE, ...) {
  # object should belong to class "lnre" if this method is called
  if (!("bootstrap" %in% names(object))) stop("no bootstrapping data available for LNRE model")
  replicates <- length(object$bootstrap)
  method <- match.arg(method)
  
  .type <- object$type
  if (missing(parm) || is.null(parm)) {
    parm <- switch(.type,
      zm=c("alpha", "B"),
      fzm=c("alpha", "B", "A"),
      gigp=c("gamma", "B", "C"))
    parm <- c(parm, "S", "X2")
  }
  
  n.parm <- length(parm)
  idx.parm <- parm %in% names(object$parm)
  idx.gof <- parm %in% names(object$gof)
  idx.S <- parm == "S"
  idx.logP <- parm == "logP"
    idx.fail <- !(idx.parm | idx.gof | idx.S | idx.logP)
  if (any(idx.fail)) stop("unknown parameter(s): ", paste(parm[idx.fail], collapse=", "))
  
  .res.list <- lapply(object$bootstrap, function (.M) {
     .res <- numeric(length=n.parm)
     if (any(idx.parm)) .res[idx.parm] <- sapply(parm[idx.parm], function (.name) .M$parm[[.name]])
     if (any(idx.gof)) .res[idx.gof] <- sapply(parm[idx.gof], function (.name) .M$gof[[.name]])
     if (any(idx.S)) .res[idx.S] <- .M$S
     if (any(idx.logP)) .res[idx.logP] <- -log10(.M$gof$p)
     names(.res) <- parm
     .res
   })
   
  .sig.level <- (1 - level) / 2  # one-sided significance level corresponding to selected confidence level
  if (method == "empirical") {
    n.outliers <- round(replicates * .sig.level) # number of "outliers" to drop from each tail
    if (n.outliers < 3) stop("insufficient bootstrap data for confidence level=", level, " (need at least ",ceiling(3 / .sig.level)," replicates)")
    .sig.level <- n.outliers / replicates # "true" significance level of the empirical interval
  }
  confint.labels <- c(sprintf("%g%%", 100 * c(.sig.level, 1 - .sig.level)), "center")

  confint.fnc <- switch(method,
    normal = function (x) {
      if (min(x) == Inf && max(x) == Inf) {
        c(Inf, Inf, Inf) # special case for infinite population diversity
      }
      else {
        .mean <- mean(x)
        .sd <- sd(x)
        .C <- -qnorm(.sig.level) # extension of confidence interval in s.d. (z-score)
        c(.mean - .C * .sd, .mean + .C * .sd, .mean)
      }
    },
    mad = function (x) {
      if (min(x) == Inf && max(x) == Inf) {
        c(Inf, Inf, Inf) # special case for infinite population diversity
      }
      else {
        .mid <- median(x)
        y <- x - .mid
        .right <- y >= 0  # data points above median
        .left <- y < 0    # data points below median
        .madR <- median(abs(y[.right])) # one-sided median absolute deviation
        .madL <- median(abs(y[.left]))
        .sdR <- .madR / qnorm(3/4) # robust MAD estimator for s.d. of Gaussian distribution
        .sdL <- .madL / qnorm(3/4) # (cf. Wikipedia article "Median absolute deviation")
        .C <- -qnorm(.sig.level) # extension of confidence interval in s.d. (z-score)
        c(.mid - .C * .sdL, .mid + .C * .sdR, .mid)
      }
    },
    empirical = function (x) {
      if (min(x) == Inf && max(x) == Inf) {
        c(Inf, Inf, Inf) # special case for infinite population diversity
      }
      else {
        y <- sort(x)
        c(y[n.outliers + 1], y[replicates - n.outliers], median(x))
      }
    }
  )
  
  .res.table <- t(apply(do.call(rbind, .res.list), 2, confint.fnc))
  rownames(.res.table) <- parm
  colnames(.res.table) <- confint.labels

  if (plot) {
    for (.p in parm) {
      x <- sapply(.res.list, function (.r) .r[.p])
      .stats <- .res.table[.p, ]
      .min <- min(x)
      .max <- max(x)
      if (.min > -Inf && .max < Inf) {
        hist(x, freq=FALSE, main=sprintf("Bootstrap of %s over %d replicates", .p, replicates), xlab=.p, xlim=c(min(.min, .stats[1]), max(.max, .stats[2])))
        abline(v=.stats[3], lwd=2, col="red")   # central value
        abline(v=.stats[1:2], lwd=1, col="red") # boundaries of confidence interval
        lines(density(x), col="blue", lwd=2)
      }
    }
  }

  .res.table
}

