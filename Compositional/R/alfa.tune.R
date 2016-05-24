################################
#### Profile log-likelihood for choosing the value of alpha
#### Fast way
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################

alfa.tune <- function(x, B = 1, ncores = 1) {
  ## x is the compositional data
  ## x must not contain any zeros
  x <- as.matrix(x)
  x <- x / rowSums(x)
  n <- nrow(x)  ## sample size
  f <- (n - 1) / n
  D <- ncol(x)  ## number of components
  d <- D - 1  ## dimensionality of the simplex
  ja <- sum( log(x) )  ## part of the Jacobian of the alpha transformation
  con <-  -n / 2 * d * log(2 * pi) - (n - 1) * d/2 + n * (d + 1/2) * log(D)

  pa <- function(a, x) {
    trans <- alfa(x, a)
    z <- trans$aff  ## the alpha-transformation
    sa <- trans$sa  ## part of the Jacobian determinant as well
    -n/2 * log( abs( det( f * cov(z) ) ) ) + (a - 1) * ja - D * sa
  }

  if (B == 1) {
    ell <- optimize(pa, c(-1, 1), x = x, maximum = TRUE )
    aff0 <- alfa(x, 0)
    z0 <- aff0$aff
    sa <- aff0$sa  ## part of the Jacobian determinant as well
    lik0 <-  -n/2 * d * log(2 * pi) - (n - 1) * d/2 -
      n/2 * log( abs( det( f * cov(z0) ) ) ) +
      n * (d + 1/2) * log(D) - ja - D * sa
    result <- c(ell$maximum, ell$objective + con, lik0)
    names(result) <- c("best alpha", "max log-lik", "log-lik at 0")

  } else {  ## bootstrap confidence intervals
    ell <- optimize(pa, c(-1, 1), x = x, maximum = TRUE )
    ab <- numeric(B)

    if (ncores == 1) {
      runtime <- proc.time()
      for (i in 1:B) {
        ind <- sample(1:n, n, replace = TRUE)
        ab[i] <- optimize(pa, c(-1, 1), x = x[ind, ], maximum = TRUE )$maximum
      }
      runtime <- proc.time() - runtime

    } else {
      runtime <- proc.time()
      cl <- makePSOCKcluster(ncores)
      registerDoParallel(cl)
      ww <- foreach::foreach( i = 1:B, .combine = rbind,
                              .export = c("alfa", "helm") ) %dopar% {
                                ind <- sample(1:n, n, replace = TRUE)
                                ab[i] <- optimize(pa, c(-1, 1), x = x[ind, ], maximum = TRUE )$maximum
                              }
      stopCluster(cl)
      ab <- as.vector( ww )
      runtime <- proc.time() - runtime
    }

    param <- c(ell$maximum, ell$objective + con, quantile( ab, c(0.025, 0.975) ) )
    names(param)[1:2] <- c("best alpha", "max log-lik")
    hist(ab, main = "Bootstrapped alpha values",
         xlab = expression( paste(alpha, " values", sep = "") ) )
    abline(v = ell$maximum, col = 3)
    abline(v = mean(ab), lty = 2, col = 4)
    message <- paste("The green is the best alpha value. The blue line is the
                     bootstrap mean value of alpha.")

    result <- list(param = param, message = message, runtime = runtime )
  }
  result
}
