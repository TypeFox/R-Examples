#' GPD Asymptotic Adjusted Information Matrix (IM) Test
#'
#' Runs the IM Test using bootstrap estimated covariance matrix. Asymptotically (in sample size) follows the F(3, bootnum - 3)
#' distribution (see reference for details).
#' @param data Data should be in vector form.
#' @param bootnum Number of bootstrap replicates for the covariance estimate.
#' @param theta Estimate for theta in the vector form (scale, shape). If NULL, uses the MLE.
#' @references Dhaene, G., & Hoorelbeke, D. (2004). The information matrix test with bootstrap-based covariance matrix estimation. Economics Letters, 82(3), 341-347.
#' @examples
#' ## Generate some data from GPD
#' x <- rgpd(200, loc = 0, scale = 1, shape = 0.2)
#' gpdImAsym(x, bootnum = 50)
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Value of theta used in the test.}
#' \item{effective_bootnum}{Effective number of bootstrap replicates used for the covariance estimate. If a
#' replicate fails to converge, it will not be used in the estimation.}
#' @export

gpdImAsym <- function(data, bootnum, theta = NULL) {
  n <- length(data)
  if(is.null(theta)) {
    fit <- tryCatch(gpdFit(data, nextremes = n, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
    if(is.null(fit))
      stop("Maximum likelihood failed to converge at initial step")
    theta <- fit$par.ests
  }
  data <- data - findthresh(data, n)
  v <- gpdImCov(data, bootnum, theta)
  eff <- v$boot_adj
  v <- v$cov
  u <- gpdInd(data, theta)
  d <- colSums(u)
  stat <- (1/n) * t(d) %*% v %*% d
  stat <- as.vector(stat)
  stat <- stat*(eff - 3) / (3*eff - 3)
  p <- 1 - pf(stat, 3, (eff - 3))
  names(theta) <- c("Scale", "Shape")
  out <- list(stat, p, theta, eff)
  names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  out
}
