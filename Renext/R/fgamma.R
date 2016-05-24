##==============================================================================
## Author: Yves Deville
##
## Find ML estimate of a two parameter gamma distribution using a sample x.
## The likelihood is concentrated with respect to the scale parameter
##==============================================================================

fgamma <- function(x,
                   check.loglik = FALSE) {

  eps <- 1e-8
  
  if (any(x <= 0)) stop("all element of x must be > 0")  
  
  n <- length(x)
  xbar <- mean(x)
  
  lx <- log(x)
  mlx <- mean(lx)

  CV <- sd(x)/xbar
  alpha <- 1/ CV / CV

  ## This is (log Lc)/n
  logLc <- function(alpha) {
    if (alpha < eps) stop("alpha must be >0")
   
    -log(gamma(alpha)) + alpha * log(alpha / xbar) +
      (alpha - 1.0) * mlx  - alpha
}
  
  res <- optimize(f = logLc, interval = c(eps, 100*alpha), maximum = TRUE)
  
  alpha.hat <- res$maximum
  beta.hat  <- xbar / alpha.hat

  loglik <- n*(-log(gamma(alpha.hat)) - alpha.hat * log(beta.hat)
               + (alpha.hat - 1.0) * mlx - alpha.hat)

  I11 <- n * trigamma(alpha.hat)
  I12 <- n / beta.hat
  I22 <- n * alpha.hat / beta.hat / beta.hat
  info <- matrix(c(I11, I12, I12, I22), nrow = 2, ncol = 2)
  covmat <- solve(info) 
  sdp <- sqrt(diag(covmat))

  estimate <- c(alpha.hat, beta.hat)
  names(estimate) <- c("shape", "scale")
  
  sdp <- sqrt(diag(covmat))
  names(sdp) <- c("shape", "scale")
  
  colnames(info) <- c("shape", "scale")
  rownames(info) <- c("shape", "scale")

  colnames(covmat) <- c("shape", "scale")
  rownames(covmat) <- c("shape", "scale")
  
  if (check.loglik) {
    check.loglik <-
      sum(stats::dgamma(x = x, shape = estimate["shape"],
                        scale = estimate["scale"],
                        log = TRUE))
  } else {
    check.loglik <- FALSE
  }
  
  list(estimate = estimate,
       sd = sdp,
       loglik = loglik,
       check.loglik = check.loglik,
       cov  = covmat,
       info = info)

}
