#' GEVr Multiplier Score Test
#'
#' Fast weighted bootstrap alternative to the parametric bootstrap procedure for the GEVr score test.
#' @param data Data should be contain n rows, each a GEVr observation.
#' @param bootnum Number of bootstrap replicates.
#' @param information To use expected (default) or observed information in the test.
#' @examples
#' x <- rgevr(500, 5, loc = 0.5, scale = 1, shape = 0.3)
#' result <- gevrMultScore(x, bootnum = 1000)
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Value of theta used in the test.}
#' @details GEVr data (in matrix x) should be of the form \eqn{x[i,1] > x[i, 2] > \cdots > x[i, r]} for each observation \eqn{i = 1, \ldots, n}.
#' @references Bader B., Yan J., & Zhang X. (2015). Automated Selection of r for the r Largest Order Statistics Approach with Adjustment for Sequential Testing. Department of Statistics, University of Connecticut.
#' @export

gevrMultScore <- function(data, bootnum, information = c("expected", "observed")) {
  data <- as.matrix(data)
  n <- nrow(data)
  R <- ncol(data)
  information <- match.arg(information)
  if(R == 1) {
    y <- tryCatch(gevrFit(data, method = "mps"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
    if(is.null(y))
      stop("MPS failed to converge at initial step")
  } else {
    y <- tryCatch(gevrFit(data[, 1:(R-1)], method = "mle"), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
    if(is.null(y))
      stop("Maximum likelihood failed to converge at initial step")
  }
  theta <- y$par.ests
  u <- gevrScore(data, theta)
  w <- colSums(u)
  if(information == "observed") {
    info <- y$varcov
  } else {
    info <- gevrFisher(data, theta)
  }
  stat <- as.vector((1/n) * t(w) %*% info %*% w)
  h <- chol(info)
  u <- u %*% t(h)
  z <- matrix(rnorm(n*bootnum, mean=0, sd=1), n, bootnum)
  z <- scale(z, center = TRUE, scale = FALSE)
  v <- t(u) %*% z
  teststat <- (1/n) * diag(t(v) %*% v)
  p <- (sum(teststat > stat) + 1) / (bootnum + 2)
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
