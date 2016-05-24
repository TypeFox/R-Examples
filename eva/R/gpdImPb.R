gpdImGen <- function(n, theta, inner) {
  x <- rgpd(n, loc = 0, scale = theta[1], shape = theta[2])
  fit1 <- tryCatch(gpdFit(x, nextremes = n, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit1)) {
    teststat <- NA
  } else {
    theta1 <- fit1$par.ests
    x <- x - findthresh(x, n)
    v <- gpdImCov(x, inner, theta1)$cov
    u <- gpdInd(x, theta1)
    d <- colSums(u)
    teststat <- (1/n) * t(d) %*% v %*% d
    teststat <- as.vector(teststat)
  }
  teststat
}


#' GPD Bootstrapped Information Matrix (IM) Test
#'
#' Runs the IM Test using a two-step iterative procedure, to boostrap the covariance estimate and critical values. See reference for details.
#' @param data Data should be in vector form.
#' @param inner Number of bootstrap replicates for the covariance estimate.
#' @param outer Number of bootstrap replicates for critical values.
#' @param allowParallel Should the outer bootstrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @references Dhaene, G., & Hoorelbeke, D. (2004). The information matrix test with bootstrap-based covariance matrix estimation. Economics Letters, 82(3), 341-347.
#' @examples
#' ## Not run
#' # x <- rgpd(200, loc = 0, scale = 1, shape = 0.2)
#' # gpdImPb(x, inner = 20, outer = 99)
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Estimate of theta for the initial dataset.}
#' \item{effective_bootnum}{Effective number of outer bootstrap replicates used (only those that converged are used).}
#' @details Warning: This test can be very slow, since the covariance estimation is nested within the outer replicates. It would be
#' recommended to use a small number of replicates for the covariance estimate (at most 50).
#' @import parallel
#' @export

gpdImPb <- function(data, inner, outer, allowParallel = FALSE, numCores = 1) {
  n <- length(data)
  fit <- tryCatch(gpdFit(data, nextremes = n, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit))
    stop("Maximum likelihood failed to converge at initial step")
  theta <- fit$par.ests
  data <- data - findthresh(data, n)
  v <- gpdImCov(data, inner, theta)$cov
  u <- gpdInd(data, theta)
  d <- colSums(u)
  stat <-  (1/n) * t(d) %*% v %*% d
  stat <- as.vector(stat)
  if(allowParallel == TRUE) {
    cl <- makeCluster(numCores)
    fun <- function(cl) {
      parSapply(cl, 1:outer, function(i,...) {gpdImGen(n, theta, inner)})
    }
    teststat <- fun(cl)
    stopCluster(cl)
  } else {
    teststat <- replicate(outer, gpdImGen(n, theta, inner))
  }
  teststat <- teststat[!is.na(teststat)]
  eff <- length(teststat)
  p <- (sum(teststat > stat) + 1) / (eff + 2)
  names(theta) <- c("Scale", "Shape")
  out <- list(stat, p, theta, eff)
  names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  out
}
