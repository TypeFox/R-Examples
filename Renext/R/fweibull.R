##==============================================================================
## Author: Yves Deville
##
## Find ML estimate of a two-parameters Weibull distribution using a sample x.
##
## The parameter 'eta' = beta^alpha is used in place of beta (= scale) where
## alpha = shape. Then 'eta' is concentrated out of the likelihood.
##==============================================================================

"fweibull" <- function(x,
                       info.observed = FALSE,
                       check.loglik = FALSE) {
  
  if (any(x <= 0)) stop("all elements in 'x' must be >0")  
  
  n <- length(x)
  xbar <- mean(x)
  
  lx <- log(x)
  mlx <- mean(lx)
  
  alpha <- 1.2825 / sd(lx)
  
  phi <- function(alpha) {
    if (alpha < 1e-8) stop("'alpha' must be > 0") 
    xa <- x^alpha
    sum(xa * lx) / sum(xa) -mlx - 1.0 / alpha
  }
  
  res <- uniroot(f = phi, interval = c(1e-8, 6 * alpha))
  
  if (abs(res$f.root) > 0.0001) {
    cat("fweibull estimation\n")
    print(res)
    stop("root not found in 'fweibull'")
  }
  
  alpha.hat <- res$root
  
  xa <- x^alpha.hat
  eta.hat <- mean(xa)
  leta <- log(eta.hat)
  
  loglik <- n * ( log(alpha.hat) - leta + (alpha.hat - 1.0) * mlx - 1.0) 
  
  if (info.observed) {
    r0 <- mean(xa)
    r1 <- mean(xa * lx)
    r2 <- mean(xa * lx * lx)

    I11 <- n * (1.0 / alpha.hat / alpha.hat + r2 / eta.hat)
    I12 <- -n * r1 / eta.hat
    I22 <- n * (2 * r0 / eta.hat -1) / eta.hat / eta.hat
    
    info <- matrix(c(I11, I12, I12, I22), nrow = 2, ncol = 2)
  } else {
    Euler <- 0.577216
    lambda1 <- 1.0 - Euler
    lambda2 <- pi * pi / 6 + Euler^2 -2 * Euler
    I11 <-  n* ( 1.0 + (leta * leta + 2 * lambda1 * leta + lambda2) ) / alpha.hat / alpha.hat 
    I12 <- -n * (lambda1 + leta) / eta.hat / alpha.hat
    I22 <- n * (1.0) / eta.hat / eta.hat
    info <- matrix(c(I11, I12, I12, I22), nrow = 2, ncol = 2)
  }

  ## derivatives of beta = scale with respect to 'alpha' and 'eta'
  dalpha <- -leta * eta.hat^(1.0 / alpha.hat) / alpha.hat/ alpha.hat
  deta <- eta.hat^(1.0 / alpha.hat - 1.0) / alpha.hat
  
  mat <- solve(info)
  sdp <- sqrt(diag(mat))

  L <- matrix(c(1, dalpha, 0, deta), nrow = 2, ncol = 2)
  cov.beta <- L%*%mat%*%t(L)

  estimate <- c(alpha.hat, eta.hat^(1/alpha.hat))
  names(estimate) <- c("shape", "scale")
  
  sdp <- sqrt(diag(cov.beta))
  names(sdp) <- c("shape", "scale")
  
  colnames(info) <- c("shape", "eta")
  rownames(info) <- c("shape", "eta")

  if (check.loglik) {
    check.loglik <-
      sum(stats::dweibull(x = x, shape = estimate["shape"],
                          scale = estimate["scale"],
                          log = TRUE))
  } else {
    check.loglik <- FALSE
  }
  
  list(estimate = estimate,
       sd = sdp,
       loglik = loglik,
       check.loglik = check.loglik,
       cov  = cov.beta,
       eta = eta.hat,
       sd.eta = sqrt(mat[2,2]),
       info = info)

}


