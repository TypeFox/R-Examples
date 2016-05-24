gevrPbGen <- function(n, R, theta, information) {
  data1 <- rgevr(n, R, theta[1], theta[2], theta[3])
  y1 <- tryCatch(gevrFit(data1, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(y1)) NA
  else gevrTestStat(y1, information)
}


#' GEVr Parametric Bootstrap Score Test
#'
#' Parametric bootstrap score test procedure to assess goodness-of-fit to the GEVr distribution.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param bootnum Number of bootstrap replicates.
#' @param information To use expected (default) or observed information in the test.
#' @param allowParallel Should the bootstrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @examples
#' ## Not run
#' ## Generate some data from GEVr
#' # x <- rgevr(200, 5, loc = 0.5, scale = 1, shape = 0.25)
#' # gevrPbScore(x, bootnum = 99)
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Initial value of theta used in the test.}
#' \item{effective_bootnum}{Effective number of bootstrap replicates (only those that converged are used).}
#' @details GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each observation \eqn{i = 1, \ldots, n}.
#' @import parallel
#' @references Bader B., Yan J., & Zhang X. (2015). Automated Selection of r for the r Largest Order Statistics Approach with Adjustment for Sequential Testing. Department of Statistics, University of Connecticut.
#' @export

gevrPbScore <- function(data, bootnum, information = c("expected", "observed"), allowParallel = FALSE, numCores = 1) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  information <- match.arg(information)
  y <- tryCatch(gevrFit(data, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(y))
    stop("Maximum likelihood failed to converge at initial step")
  theta <- y$par.ests
  stat <- gevrTestStat(y, information)
  if(allowParallel == TRUE) {
    cl <- makeCluster(numCores)
    fun <- function(cl) {
      parSapply(cl, 1:bootnum, function(i,...) {gevrPbGen(n, R, theta, information)})
    }
    teststat <- fun(cl)
    stopCluster(cl)
  } else {
    teststat <- replicate(bootnum, gevrPbGen(n, R, theta, information))
  }
  teststat <- teststat[!is.na(teststat)]
  eff <- length(teststat)
  p <- (sum(teststat > stat) + 1) / (eff + 2)
  out <- list(stat, p, theta, eff)
  names(out) <- c("statistic", "p.value", "theta", "effective_bootnum")
  out
}
