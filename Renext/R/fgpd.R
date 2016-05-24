##**********************************************************************
## Maximum Likelihood estimation of a Generalised Pareto Distibution
## by likelihood concentration.
##
## TODO : situations with one known parameter
##**********************************************************************

fGPD <- function(x,
                 info.observed = FALSE,
                 shapeMin = -0.8,
                 dCV = 1e-4,
                 cov = TRUE,
                 trace = 0) {
  
  parnames <- c("shape", "scale")
  n <- length(x)
  M1 <- mean(x)
  CV <- sqrt(1 - 1/n) * sd(x) / M1

  if (any(x <= 0)) stop("all elements in 'x' must be > 0")
  if (shapeMin >= 0) stop("'shapeMin' must be negative")
  if (shapeMin < -1) warning("'shapeMin < -1': non-identifiable model")
  if (dCV >= 1e-2) warning("'dCV' should be small < 0.01")
  
  if (CV > 1.0 + dCV) {
      ## Lomax case
      res.lomax <- flomax(x,
                          info.observed = info.observed, cov = cov)
      res.gpd <- lomax2gpd(parLomax = res.lomax$estimate, vcovLomax = res.lomax$cov)
      loglik <- res.lomax$loglik
      if (!is.null(vcov <- attr(res.gpd, "vcov"))) {
          dloglik <- attr(res.gpd, "jacobian")[ , c("shape", "scale")] %*% res.lomax$dloglik
          sd.gpd <- sqrt(diag(vcov))
      } else {
          dloglik <- NULL
          sd.gpd <- NULL
      }
  } else if (CV < 1.0 - dCV) {
      ## Maxlo case
      res.maxlo <- fmaxlo(x, shapeMin = - 1 / shapeMin,
                          info.observed = info.observed, cov = cov)
      res.gpd <- maxlo2gpd(parMaxlo = res.maxlo$estimate,
                           vcovMaxlo = res.maxlo$cov)
      loglik <- res.maxlo$loglik
      if (!is.null(vcov <- attr(res.gpd, "vcov"))) {
          dloglik <- attr(res.gpd, "jacobian")[ , c("shape", "scale")] %*% res.maxlo$dloglik
          sd.gpd <- sqrt(diag(vcov))
      } else {
          return(list(estimate = res.gpd[1L:2L], CV = CV, loglik = loglik, cvg = TRUE))
      }
  } else {
      ## exponential case. Note that the variance of the estimated scale
      ## is greater than that of the one parameter exponential distribution
      if (trace) {
         cat(sprintf("CV = %6.4f close to 1: exponential\n", CV))
      }
      res.gpd <- c("scale" = M1, "shape" = 0.0)
      loglik <- -n * (1 + log(M1)) 
      if (cov) {
          dloglik <- c("scale" = 0.0, "shape" = 0.0)
          vcov <- matrix(c(2 * M1^2, -M1, -M1, 1) / n, nrow = 2, ncol = 2)
          colnames(vcov) <- rownames(vcov) <- c("shape", "scale")
          rn <- sqrt(n)
          sd.gpd <- c("scale" = M1 * sqrt(2) / rn,  shape = 1 / rn)
      } else {
          return(list(estimate = res.gpd[1L:2L], CV = CV, loglik = loglik, cvg = TRUE))
      }

  }
  
  list(estimate = res.gpd[1L:2L],
       CV = CV,
       loglik = loglik,
       dloglik = dloglik,
       sd = sd.gpd,
       cov = vcov,
       cvg = TRUE)
    
}

##==============================================================================
## Author: Yves Deville
##
## Find ML estimate of a one parameter GPD (shape 'xi' fixed)
##
## The lilelihood is easily maximised in that case
##==============================================================================

fgpd1 <- function(x, shape = 0.2, plot = FALSE) {

  xi <- shape
  names(xi) <- NULL
  if (xi <= -1.0) stop("parameter 'shape' must be > -1")
  
  n <- length(x)
  
  if (any(is.na(x)) || any(x < 0) ) stop("'x' elements must be > 0 and non NA") 
  
  ## exponential case...
  if (xi == 0.0) {
    sigma.hat <- mean(x)
    cov0 <- sigma.hat^2 / n
    return(list(estimate = c(scale = sigma.hat),
                sd = sigma.hat / sqrt(n),
                loglik = -log(sigma.hat) - n,
                cov  = cov0,
                info = 1 / cov0))
  }
  
  ## not used for now
  dlogL1 <- function (sigma) {
    xmod <- x / sigma
    ( -n + (xi + 1.0) * sum(xmod / (1.0 + xi * xmod) ) ) / sigma 
  }
  
  logL1 <- function (sigma) {
    xmod <- x / sigma
    -n * log(sigma)  - (xi + 1) * sum( log(1 + xi * xmod) ) / xi
  }

  if (xi < 0 ) {
    ## when  xi < 0, mind the support
    interv <- c(-xi * max(x), max(x))
  } else {
    ## when  xi > 0, the logL is increasing at min(x)
    ## and decreasing at max(x)
    interv <- range(x)
  }

  res <- optimize(f = logL1, interval = interv, maximum = TRUE)
  sigma.hat <- res$maximum
  
  loglik <- res$objective

  if (xi > -0.5) {
    info <- n / (1 + 2 * xi) / sigma.hat / sigma.hat
    cov0 <- 1 / info
    sdp <- sqrt(cov0)
  } else {
    warning("'shape' is <= -0.5. No variance provided for the estimator of 'shape'") 
    info <- NULL
    cov0 <- NULL
    sdp <- NULL
  }
     
  if (plot) {
    sigmas <- seq(from = interv[1], to = max(x), length.out = 500)
    lL <- sapply(sigmas, logL1)
    plot(x = sigmas, y = lL, type = "l", cex = 0.8,
         main = sprintf("xi = %4.2f", xi))
    abline(v = res$maximum, col = "red")
    abline(v = res$maximum + c(-1, 1) * sdp, col = "red", lty = "dashed")
    ## abline(v = sigma, col = "SpringGreen3")
  }
  
  list(estimate = c(scale = sigma.hat),
       sd = sdp,
       loglik = loglik,
       cov  = cov0,
       info = info)

}
  
