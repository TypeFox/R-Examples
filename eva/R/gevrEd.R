#' GEVr Entropy Difference Test
#'
#' Goodness-of-fit test for GEVr using the difference in likelihood between GEVr and GEV(r-1).
#' This can be used sequentially to test for the choice of r.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param theta Estimate for theta in the vector form (loc, scale, shape). If NULL, uses the MLE from the full data.
#' @examples
#' ## This will test if the GEV2 distribution fits the data.
#' x <- rgevr(100, 2, loc = 0.5, scale = 1, shape = 0.5)
#' result <- gevrEd(x)
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Estimate of theta using the top r order statistics.}
#' @details GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each
#' observation \eqn{i = 1, \ldots, n}. The test uses an asymptotic normality result based on the expected
#' entropy between the GEVr and GEV(r-1) likelihoods. See reference for detailed information. This test can be
#' used to sequentially test for the choice of r, implemented in the function `gevrSeqTests'.
#'
#' @references Bader B., Yan J., & Zhang X. (2015). Automated Selection of r for the r Largest Order Statistics Approach with Adjustment for Sequential Testing. Department of Statistics, University of Connecticut.
#' @import stats graphics
#' @export
gevrEd <- function(data, theta = NULL) {
  data <- as.matrix(data)
  R <- ncol(data)
  if(R == 1)
    stop("R must be at least two")
  n <- nrow(data)
  if(is.null(theta)) {
    y <- tryCatch(gevrFit(data, method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
    if(is.null(y))
      stop("Maximum likelihood failed to converge at initial step")
    theta <- y$par.ests
  }
  Diff <- dgevr(data[, 1:R], loc = theta[1], scale = theta[2], shape = theta[3], log.d = TRUE) -
          dgevr(data[, 1:(R-1)], loc = theta[1], scale = theta[2], shape = theta[3], log.d = TRUE)
  EstVar <- sum((Diff - mean(Diff))^2) / (n-1)
  FirstMom  <- -log(theta[2]) - 1 + (1 + theta[3]) * digamma(R)
  Diff <- sum(Diff) / n
  Diff <- sqrt(n)*(Diff  - FirstMom) / sqrt(EstVar)
  p.value <- 2*(1-pnorm(abs(Diff)))
  out <- list(as.numeric(Diff), as.numeric(p.value), theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
