#' @rdname LambertW-utils
#' @description
#' \code{qqLambertW} computes and plots the sample quantiles of the data
#' \code{y} versus the theoretical Lambert W \eqn{\times} \eqn{F} theoretical
#' quantiles given \eqn{\theta}.
#' @export
qqLambertW <- function(y, distname, theta = NULL, beta = NULL, gamma = 0, delta = 0, alpha = 1, 
                       plot.it = TRUE, use.mean.variance = TRUE, ...) {
  
  
  stopifnot(is.numeric(y),
            length(y[!is.na(y)]) > 0)
  
  check_distname(distname)
  if (is.null(theta)) {
    warning("Please specify parameters by passing a list",
            "to the 'theta' argument directly.\n",
            "Specifying parameters by alpha, beta, gamma, delta will be",
            "deprecated.")
    theta <- list(beta = beta, alpha = alpha, gamma = gamma, delta = delta)
  } 
  theta <- complete_theta(theta)
  check_theta(theta = theta, distname = distname)
  
  xlab <- "Theoretical Quantiles"
  ylab <- "Sample Quantiles"
  
  main <- paste("Lambert W x", distname, "QQ plot")
  
  y <- y[!is.na(y)]
  nn <- length(y)
  
  p.n <- ppoints(nn)
  x <- qLambertW(p.n, theta = theta, distname = distname,
                 use.mean.variance = use.mean.variance)
  sorted.x <- sort(x)
  sorted.y <- sort(y)
  if (plot.it) {
    plot(sorted.x, sorted.y, main = main, xlab = xlab, ylab = ylab, ...)
    abline(0, 1)
  }
  invisible(list(x = sorted.x, y = sorted.y))
} 