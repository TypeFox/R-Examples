# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

resample <- function(data, resampleFun, sampler, R = 10000,
                     seed = NULL,
                     statisticNames = NULL,
                     block.size = 100, trace = FALSE, ...,
                     observedIndices = 1:n,
                     call = match.call())
{
  # Perform resampling. This is the workhorse function, called by bootstrap
  # and other resampling functions to actually do the resampling.
  #
  # Args:
  #   data:            vector, matrix, or data frame
  #   resampleFun:     function(data, ii), created by .resampleMakeFun
  #   sampler:         function(n, R, ...) to generate samples
  #   R:               number of replications
  #   seed:            old value of .Random.seed, or argument to set.seed
  #   statisticNames:  names used for printing, character vector of length 'd'
  #   block.size:      replicates are done 'block.size' at a time
  #   trace:           logical, if TRUE an indication of progress is printed
  #   ...:             additional arguments passed to sampler
  #   observedIndices  vector of indices that correspond to the observed value.
  #                    This is used by permutationTest2
  #   call:            The call to the function that is calling resample.

  .resampleSetSeed(seed)
  oldSeed <- .Random.seed

  dimData <- dim(data)
  n <- IfElse(is.null(dimData), length(data), dimData[1])

  # Observed value
  observed <- resampleFun(data, observedIndices)
  p <- length(observed)
  if(p) { # p = 0 for graphical bootstrapping
    if(!is.null(dim(observed)))
      dim(observed) <- NULL
    if(!is.null(statisticNames) && length(statisticNames) != p)
      stop("statisticNames has the wrong length, is ",
           length(statisticNames),
           ", but observed has length ", p)
    ## Preference for names:
    # supplied names
    # names(observed)
    # names constructed from statistic
    # stat1 stat2 ...
    # TODO: construct names that combine data & statistic. E.g. for colMeans.
    # TODO: when passing args.stat, reflect that in names
    if(is.null(statisticNames))
      statisticNames <- names(observed)
    if(is.null(statisticNames)) {
      temp <- call$statistic
      temp2 <- IfElse(is.call(temp) && temp[[1]] == "function", "stat",
                      # e.g. c(mean = mean(x), sd = sd(x))
                      is.call(temp) && p > 1, "stat",
                      # e.g. mean or mean(X)
                      deparse(temp))
      statisticNames <- paste0(temp2, IfElse(p == 1, NULL, 1:p))
    }
    #if(is.null(statisticNames))
    #  statisticNames <- paste0("stat", 1:p)
    temp <- which(statisticNames == "")
    if(length(temp))
      statisticNames[temp] <- paste0("stat", temp)
    names(observed) <- statisticNames

    replicates <- matrix(NA * vector(storage.mode(observed), R * p),
                         nrow = R, ncol = p,
                         dimnames = list(NULL, statisticNames))
  }
  failures <- NULL

  # Do calculations, mostly block.size at a time.
  # Generate indices using sampler, and statistics using resampleFun.
  nBlocks <- ceiling(R / block.size)
  B1 <- 0
  for(iBlock in 1:nBlocks) {
    B0 <- B1 # Number of replications already performed
    B1 <- min(R, iBlock * block.size) # Last replication to do in this block
    if(trace)
      cat("Calculating replications ", B0 + 1, ":", B1, "\n", sep = "")
    BB <- B1 - B0
    indices <- sampler(n, BB, ...)
    for(j in 1:BB) {
      theta <- try(resampleFun(data, indices[, j]), silent = TRUE)
      if(is(theta, "try-error")) {
        failures <- c(failures, B0 + j)
      } else if(p) {
        replicates[B0 + j, ] <- theta
      }
    }
  }
  if(!p)
    replicates <- NULL
  result <- list(observed = observed,
                 replicates = replicates,
                 n = n, p = p, R = R,
                 seed = oldSeed, call = call)
  result$failures <- failures
  if(!missing(observedIndices)) {
    # Working with a subset of a larger dataset.
    result$n <- length(observedIndices)
    result$nCombined <- n
  }
  class(result) <- "resample"
  return(result)
}


print.resample <- function(x, ...) {
  # This is mostly used for objects with $call and $stats, like "bootstrap"
  if(!is.null(x$call))
    cat0n("Call:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"))
  catn("Replications:", x$R)
  if(is.null(x$ratio)) { # one-sample application
    catn("\nSummary Statistics:")
  } else {
    catn("Two samples, sample sizes are", x$n[1], x$n[2])
    catn("\nSummary Statistics for the",
         IfElse(x$ratio, "ratio", "difference"),
         "between samples 1 and 2:")
  }
  print(IfElse(is.null(x$stats),
               .BootStats(x)[1:3], # Observed SE Mean, not Bias
               x$stats), ...)
  if(length(x$failures))
    catn(length(x$failures), "replicates failed")
  invisible(x)
}

hist.resample <- function(x, ..., resampleColumns = 1:x$p, xlim = NULL,
                          xlab = NULL, main = "", col = "blue", border = 0,
                          breaks = "FD", showObserved = TRUE,
                          legend = TRUE, args.legend = NULL) {
  statNames <- names(x$observed)
  nullXlim <- is.null(xlim)
  for(j in resampleColumns) {
    xj <- x$replicates[, j]
    if(nullXlim)
      xlim <- range(xj, x$observed[j], na.rm = TRUE)
    hist(xj, probability = TRUE,
         xlab = IfElse(is.null(xlab), statNames[j], xlab),
         xlim = xlim, ..., main = main, col = col, border = border,
         breaks = breaks)
    den <- density(xj, from = xlim[1], to = xlim[2])
    lines(den)
    if(showObserved) {
      om <- c(x$observed[j], x$stats$Mean[j]) # observed & mean
      points(om, rep(0, 2), pch = 1:2)
      tops <- approx(den, xout = om)$y
      segments(om, y0 = 0, y1 = pmax(tops, .2 * max(den$y)), lty = 1:2)
      if(legend) {
        args.legend <- c(args.legend,
                         list(x = "topright", legend = c("Observed", "Mean"),
                              lty = 1:2, pch = 1:2, box.col = 0))
        args.legend <- args.legend[!duplicated(names(args.legend))]
        do.call("legend", args.legend)
      }
    }
  }
  invisible(NULL)
}

plot.resample <- function(x, ...) hist.resample(x, ...)

qqnorm.resample <- function(y, ..., resampleColumns = 1:y$p,
                            ylab = NULL, pch = if(y$R < 100) 1 else ".") {
  statNames <- names(y$observed)
  for(j in resampleColumns) {
    yj <- y$replicates[, j]
    qqnorm(yj, ..., pch = pch,
           ylab = IfElse(is.null(ylab), statNames[j], ylab))
  }
  invisible(NULL)
}

quantile.resample <- function(x, ...) {
  # Use type = 6 for wider quantiles, better accuracy
  # Return 1 row per statistic, and one column for each prob
  t(apply(x$replicates, 2, Quantile, ...))
}

Quantile <- function(x, ..., type = 6){
  # Use type = 6 for better accuracy.
  quantile(x, ..., type = type)
}

if(FALSE) { # manual testing code
  r <- bootstrap((xDF), colMeans) # xDF defined in bootstrap.R
  plot(r, resampleColumns = 2)
  qqnorm(r, resampleColumns = 2)
  quantile(r, probs = c(.025, .975))
  limits.percentile(r, probs = c(.025, .975))

  r <- bootstrap2(treatment = t9, x9, mean)
  plot(r)

  r <- permutationTest2(treatment = t9, x9, mean)
  plot(r)
}
