#' Sequential Tests for the GEVr Model
#'
#' Sequentially performs the entropy difference (ED) test or the multiplier or parametric bootstrap score tests for the GEVr model.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param method Which test to run: ED test (ed), multiplier (multscore) or parametric bootstrap (pbscore) score test.
#' @param bootnum If method equals 'pbscore' or 'multscore', the number of bootstrap simulations to use.
#' @param information To use expected (default) or observed information in the score tests.
#' @param allowParallel If method equals 'pbscore', should the parametric boostrap procedure be run in parallel or not. Defaults to false.
#' @param numCores If allowParallel is true, specify the number of cores to use.
#' @examples
#' x <- rgevr(200, 5, loc = 0.5, scale = 1, shape = 0.25)
#' gevrSeqTests(x, method = "ed")
#' @return
#' Function returns a dataframe containing the test statistics, estimates, and p-value results of the sequential tests.
#'
#' \item{r}{Value of r to be tested.}
#' \item{p.values}{Raw p-values from the individual tests at each value of r.}
#' \item{ForwardStop}{Transformed p-values according to the ForwardStop stopping rule.}
#' \item{StrongStop}{Transformed p-values according to the StrongStop stopping rule.}
#' \item{statistic}{Returned test statistics of each individual test.}
#' \item{est.loc}{Estimated location parameter for the given r.}
#' \item{est.scale}{Estimated scale parameter for the given r.}
#' \item{est.shape}{Estimated shape parameter for the given r.}
#' @details GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each observation \eqn{i = 1, \ldots, n}.
#' See function `pSeqStop' for details on transformed p-values.
#' @export

gevrSeqTests <- function(data, bootnum = NULL, method = c("ed", "pbscore", "multscore"), information = c("expected", "observed"),
                              allowParallel = FALSE, numCores = 1) {
  data <- as.matrix(data)
  R <- ncol(data)
  method <- match.arg(method)
  if(method != "ed") {
    if(is.null(bootnum))
      stop("Must enter the number of bootstrap replicates!")
    information <- match.arg(information)
    result <- matrix(0, R, 8)
    for(i in 1:R) {
      result[i, 1] <- i
      if(method == "multscore")
        fit <- gevrMultScore(data[, 1:i], bootnum, information)
      if(method == "pbscore")
        fit <- gevrPbScore(data[, 1:i], bootnum, information, allowParallel, numCores)
      result[i, 2] <- fit$p.value
      result[i, 5] <- fit$statistic
      result[i, 6:8] <- fit$theta
    }
  } else {
    if(R == 1)
      stop("R must be at least two")
    result <- matrix(0, R-1, 8)
    for(i in 2:R) {
      result[i-1, 1] <- i
      fit <- gevrEd(data[, 1:i])
      result[i-1, 2] <- fit$p.value
      result[i-1, 5] <- fit$statistic
      result[i-1, 6:8] <- fit$theta
    }
  }
  result[, 3] <- rev(pSeqStop(rev(result[, 2]))$ForwardStop)
  result[, 4] <- rev(pSeqStop(rev(result[, 2]))$StrongStop)
  colnames(result) <- c("r", "p.values", "ForwardStop", "StrongStop", "statistic", "est.loc", "est.scale", "est.shape")
  as.data.frame(result)
}
