################################
#### Multivariate t distribution
#### Tsagris Michail 10/2012
#### mtsagris@yahoo.gr
################################

multivt <- function(y, plot = FALSE) {
  ## the next mvt function is for the appropriate
  ## degrees of freedom
  ## y contains the data
  y <- as.matrix(y)  ## makes sure y is a matrix
   mvt <- function(y, v) {
    ## the next function 'a' estimates the mean and covariance for given
    ## degeees of freedom. It's a built-in function
    a <- MASS::cov.trob(y, nu = v)
    se <- a$cov
    n <- nrow(y)
    p <- ncol(y)
    me <- as.vector(a$center)
    f <- n * lgamma( (v + p)/2 ) - n * lgamma(v/2) - 0.5 * n * p *
    log(pi * v) - 0.5 * n * log( det(se) ) - 0.5 * (v + p) *
    sum( log(1 + mahalanobis(y, me, se)/v) )
    f
   }
  b <- optimize(mvt, c(0.9, 20000), y = y, maximum = T)
  dof <- b$maximum
  loglik <- b$objective
  ## df is the optimal degrees of freedom
  ## if the df is a big number, then a multivariate normal is fine as well
  result <- cov.trob(y, nu = dof)  ## the center and covariance matrix
  ## will be calculated based on the optimal degrees of freedom
  ## the classical mean and covariance are given in the results
  ## for comparison pruposes
  apotelesma <- list(center = result$center, scatter = result$cov,
  df = dof, loglik = loglik, mesos = colMeans(y), covariance = cov(y))
  if (plot == TRUE) {
    lik <- deg <- seq(max(1, dof - 20), dof + 20, by = 0.1)
    for ( i in 1:length(deg) )  lik[i] <- mvt(y, deg[i])
    plot(deg, lik, type = "l", xlab = "Degrees of freedom",
    ylab = "Log likelihood")
    b <- max(lik) - 1.92
    abline(h = b, col = 2)
    a1 <- min(deg[lik >= b])
    a2 <- max(deg[lik >= b])
    abline(v = a1, col = 3, lty = 2)
    abline(v = a2, col = 3, lty = 2)
    conf <- c(a1, a2)
    names(conf) <- c("2.5%", "97.5%")
    apotelesma <- list(center = result$center, scatter = result$cov,
    df = dof, conf = conf, loglik = loglik, mesos = colMeans(y),
    covariance = cov(y))
  }
  apotelesma
}
