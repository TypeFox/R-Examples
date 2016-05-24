#' GPD Multiplier Score Test
#'
#' Fast weighted bootstrap alternative to the parametric bootstrap procedure for the Generalized Pareto score test.
#' @param data Data should be in vector form.
#' @param bootnum Number of bootstrap replicates.
#' @param information To use expected (default) or observed information in the test.
#' @examples
#' x <- rgpd(100, loc = 0, scale = 1, shape = 0.25)
#' gpdMultScore(x, bootnum = 1000)
#' @return
#' \item{statistic}{Test statistic.}
#' \item{p.value}{P-value for the test.}
#' \item{theta}{Value of theta used in the test.}
#' @export

gpdMultScore <- function(data, bootnum, information = c("expected", "observed")) {
  n <- length(data)
  information <- match.arg(information)
  y <- tryCatch(gpdFit(data, method = "mps", nextremes = n), error = function(w) {return(NULL)}, warning = function(w) {return(NULL)})
  if(is.null(y))
    stop("MPS failed to converge at initial step")
  theta <- y$par.ests
  data <- data - findthresh(data, n)
  u <- gpdScore(data, theta)
  w <- colSums(u)
  if(information == "observed") {
    info <- y$varcov
  } else {
    info <- gpdFisher(n, theta)
  }
  stat <- t(w) %*% info %*% w
  stat <- as.vector(stat)
  h <- chol(info)
  u <- u %*% t(h)
  z <- matrix(rnorm(n*bootnum, mean = 0, sd = 1), n, bootnum)
  z <- scale(z, center = TRUE, scale = FALSE)
  v <- t(u) %*% z
  teststat <- diag(t(v) %*% v)
  p <- (sum(teststat > stat) + 1) / (bootnum + 2)
  names(theta) <- c("Scale", "Shape")
  out <- list(stat, p, theta)
  names(out) <- c("statistic", "p.value", "theta")
  out
}
