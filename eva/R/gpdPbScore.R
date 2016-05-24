gpdPbGen <- function(n, theta, information) {
  data1 <- rgpd(n, loc = 0, scale = theta[1], shape = theta[2])
  fit1 <- tryCatch(gpdFit(data1, nextremes = n, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit1)) NA
  else gpdTestStat(fit1, information)
}


#' GPD Parametric Bootstrap Score Test
#'
#' Parametric bootstrap score test procedure to assess goodness-of-fit to the Generalized Pareto distribution.
#' @param data Data should be in vector form.
#' @param bootnum Number of bootstrap replicates.
#' @param information To use expected (default) or observed information in the test.
#' @param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @examples
#' ## Generate some data from GPD
#' x <- rgpd(200, loc = 0, scale = 1, shape = 0.2)
#' gpdPbScore(x, bootnum = 100)
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Estimated value of theta for the initial data.}
#' \item{effective_bootnum}{Effective number of bootstrap replicates (only those that converged are used).}
#' @import parallel
#' @export

gpdPbScore <- function(data, bootnum, information = c("expected", "observed"), allowParallel = FALSE, numCores = 1) {
  n <- length(data)
  information <-  match.arg(information)
  fit <- tryCatch(gpdFit(data, nextremes = n, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(fit))
    stop("Maximum likelihood failed to converge at initial step")
  theta <- fit$par.ests
  stat <- gpdTestStat(fit, information)
  if(allowParallel == TRUE) {
    cl <- makeCluster(numCores)
    fun <- function(cl) {
      parSapply(cl, 1:bootnum, function(i,...) {gpdPbGen(n, theta, information)})
    }
    teststat <- fun(cl)
    stopCluster(cl)
  } else {
    teststat <- replicate(bootnum, gpdPbGen(n, theta, information))
  }
  teststat <- teststat[!is.na(teststat)]
  eff <- length(teststat)
  p <- (sum(teststat > stat) + 1) / (eff + 2)
  names(theta) <- c("Scale", "Shape")
  out <- list(stat, p, theta, eff)
  names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  out
}
