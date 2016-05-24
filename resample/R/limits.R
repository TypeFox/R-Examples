# Copyright 2014 Google Inc. All rights reserved.
#
# Use of this source code is governed by a BSD-style
# license that can be found in the LICENSE file or at
# http://opensource.org/licenses/BSD-3-Clause

CI.percentile <- function(x, probs = c(.025, .975),
                              expand = TRUE, ...) {
  # Bootstrap percentile confidence interval
  # Args:
  #   x: a bootstrap or bootstrap2 object
  #   probs: confidence levels
  #   expand: logical, if TRUE then use modified percentiles for better
  #              small-sample accuracy.
  probs2 <- IfElse(expand, ExpandProbs(probs, min(x$n)), probs)
  result <- Quantile(x, probs = probs2, ...)
  dimnames(result) <- list(names(x$observed), .FormatProbs(probs))
  return(result)
}
# TODO: separate calibration for t and SE effects

limits.percentile <- function(...) {
  warning("limits.percentile is deprecated, use CI.percentile instead")
  CI.percentile(...)
}

ExpandProbs <- function(probs, n) {
  # Modify the probs to adjust for two factors:
  # Bootstrap distributions are too narrow by a factor of sqrt((n-1)/n)
  # Percentile limits correspond to estimate +- z_alpha instead of t_alpha s.
  # Find probs2 such that qnorm(probs2) * sqrt((n-1)/n) = qt(probs, n-1)
  pnorm(qt(probs, n-1) * sqrt(n / (n-1)))
}

.FormatProbs <- function(probs) {
  # Format quantiles, e.g. 0.025 becomes 2.5%
  paste0(formatC(100 * probs, format = "fg", width = 1,
                 max(2L, getOption("digits"))), "%")
}

CI.t <- function(x, probs = c(.025, .975)) {
  # T interval with bootstrap SE
  obs <- x$observed
  SE <- x$stats$SE
  df <- min(x$n) - 1
  result <- rep(obs, each = length(probs)) + qt(probs, df = df) %o% SE
  dimnames(result) <- list(.FormatProbs(probs), names(obs))
  t(result)
}

limits.t <- function(...) {
  warning("limits.t is deprecated, use CI.t instead")
  CI.t(...)
}


CI.bootstrapT <- function(x, probs = c(.025, .975)) {
  # Bootstrap t interval
  # This assumes that the first dimension of the statistic
  # is an estimate, and the second is proportional to a SE for the estimate.
  # E.g. for bootstrapping the mean, they could be the mean and s.
  # Args:
  #   x: bootstrap object
  tAlpha <- Quantile((x$replicates[, 1] - x$observed[1]) / x$replicates[, 2],
                     1-probs)
  result <- x$observed[1] - tAlpha * x$observed[2]
  matrix(result, 1, dimnames = list(names(x$observed[1]), .FormatProbs(probs)))
}
# S+Resample has bootstrapT, with a pivot argument.

limits.bootstrapT <- function(...) {
  warning("limits.bootstrapT is deprecated, use CI.bootstrapT instead")
  CI.bootstrapT(...)
}


CI.bca <- function(x, probs = c(.025, .975),
                       expand = TRUE, L = NULL, ...) {
  # Bootstrap BCa confidence interval
  # Args:
  #   x: a bootstrap object
  #   probs: confidence levels
  #   L:     empirical influence function
  if(is(x, "bootstrap2"))
    stop("bootstrap2 objects are not currently supported")

  skewness <- function(x) {
    x <- x - mean(x)
    mean(x^3) / mean(x^2)^1.5
  }
  if(is.null(L)){
    Call <- x$call
    Call[[1]] <- as.name("jackknife")
    Call$R <- NULL
    Call$seed <- NULL
    Call$sampler <- NULL
    Call$block.size <- NULL
    jackObject <- eval(Call, sys.parent())
    L <- -jackObject$replicates
  }
  a <- apply(L, 2, skewness) / (6 * sqrt(nrow(L)))
  w <- qnorm(colMeans(x$replicates < rep(x$observed, each = x$R)))
  probs2 <- IfElse(expand, ExpandProbs(probs, min(x$n)), probs)
  zalpha <- qnorm(probs2)
  # For now probs in rows, statistics in columns
  zalpha <- matrix(zalpha, nrow = length(probs), ncol = x$p)
  probs3 <- pnorm(w + (w + zalpha) / (1 - a * (w + zalpha)))
  result <- probs3 * NA
  for(j in 1:x$p) {
    result[, j] <- Quantile(x$replicates[, j], probs3[, j])
  }
  dimnames(result) <- list(.FormatProbs(probs), names(x$observed))
  t(result)
}
# TODO: test that.
# TODO: Add expand argument.
# TODO: Support two-sample applications.



if(FALSE) {
  x9 <- 1:9
  xDF <- data.frame(X = x9, Y = 2*x9)

  r <- bootstrap(x9, mean, R=100)
  CI.t(r)
  CI.percentile(r)
  CI.bca(r)

  CI.t(r, probs = .9)
  CI.percentile(r, probs = .9)
  CI.bca(r, probs = .9)


  r <- bootstrap(xDF, colMeans, R=100)
  CI.t(r)
  CI.percentile(r)
  CI.bca(r)

  r <- bootstrap2(treatment = t9, xDF, colMeans, R=100)
  CI.t(r)
  CI.percentile(r)

  r <- bootstrap((1:10)^2, c(mean(data), sd(data)), R=100)
  CI.bootstrapT(r)
  CI.bootstrapT(r, probs = .9)
}


if(FALSE){ # Development work, using examples in bootstrap/howToTeach
  for(i in list.files(path = "~/resample/resample/R", pattern = "*.R$",
                      full.names = TRUE)) source(i)

  jackBasic <- jackknife(Basic, mean)
  jack2 <- jackknife(Basic, function(x) mean((x-mean(x))^2))
  jack2
  # Compare that bias to
  sd(Basic)^2 - mean((Basic-mean(Basic))^2)
  CI.t(jackBasic)
  CI.t(jack2)

  boot <- bootstrap(Basic, c(mean = mean(Basic), median = median(Basic)))
  jack <- jackknife(Basic, c(mean = mean(Basic), median = median(Basic)))
  CI.percentile(bootBasic)
  CI.percentile(boot)
  CI.percentile(bootBasic, expand = TRUE)

  CI.bca(bootBasic)
  CI.bca(boot)
}
