################################
#### Wood's bimodal distribution on the sphere
#### Tsagris Michail 1/2016 
#### mtsagris@yahoo.gr
#### References: Andrew T.A. Wood (1982), JRSSC, 31(1): 52-58
#### A bimodal distribution on the sphere 
################################

wood.mle <- function(y) {
  ## y is a two column matrix, where the first column is the latitude and 
  ## the second is the longitude, all expressed in degrees
  y <- as.matrix(y)
  y[, 1] <- 90 - y[, 1] ## we want the co-latitude
  y <- y / 180 * pi
  x <- cbind( sin(y[, 1]) * cos(y[, 2]), sin(y[, 1]) * sin(y[, 2]), cos(y[, 1]) )
  mle <- function(param) {
    gam <- param[1]  ;  del <- param[2]
    m1 <- c( cos(gam) * cos(del), cos(gam) * sin(del), - sin(gam) )
    m2 <- c( - sin(del), cos(del), 0 )
    m3 <- c( sin(gam) * cos(del), sin(gam) * sin(del), cos(gam) )
    a1 <- as.vector( x %*% m1 )
    a2 <- as.vector( x %*% m2 )
    a3 <- as.vector( x %*% m3 )
    u <- sum(a3)
    v <- sum( ( a1^2 - a2^2 ) / sqrt( 1 - a3^2 ) )  
    w <- 2 * sum( a1 * a2  / sqrt( 1 - a3^2 ) )
    f <-  - ( u^2 + v^2 + w^2 ) 
    f
  }  
  lik <- function(k) {
    - n * log(2 * pi) - n * log( ( exp(k) - exp(-k) ) / k ) + 
    k * ( u * cos(a) + ( v * cos(b) + w * sin(b) ) * sin(a) )
  }
  ini <- colMeans(y)
  mod <- optim( ini, mle )
  mod <- optim(mod$par, mle, hessian = TRUE) 
  gam <- mod$par[1]  ;  del <- mod$par[2]
  m1 <- c( cos(gam) * cos(del), cos(gam) * sin(del), - sin(gam) )
  m2 <- c( - sin(del), cos(del), 0 )
  m3 <- c( sin(gam) * cos(del), sin(gam) * sin(del), cos(gam) )
  a1 <- as.vector( x %*% m1 )
  a2 <- as.vector( x %*% m2 )
  a3 <- as.vector( x %*% m3 )
  u <- sum(a3)
  v <- sum( ( a1^2 - a2^2 ) / sqrt( 1 - a3^2 ) )  
  w <- 2 * sum( a1 * a2  / sqrt( 1 - a3^2 ) )
  a <- atan2( sqrt(v^2 + w^2), u ) 
  b <- atan2( w, v )
  n <- nrow(y)
  ka <- optimize( lik, c(0, 1000), maximum = TRUE )
  k <- ka$maximum
  info <- c( c( mod$par, b) / pi * 180, ka$objective)
  names(info) <- c("gamma", "delta", "beta", "log-lik")
  se.mod <- sqrt( diag( solve(mod$hessian) ) )
  se.a <- sqrt( (1 - exp( - 2 * k)) / (1 + exp(- 2* k ) ) *( 1 / (n * k) ) )
  se.b <- 1 / ( k * ( 1 + exp(-2 * k) ) / ( 1 - exp(- 2 * k) ) * sin(a)^2 )
  se.k <- sqrt( 1 / ( 1 / k^2 - 2 / ( exp(k) - exp(- k) ) ) ) / sqrt(n)
  confgam <- c(gam - qnorm(0.975) * se.mod[1], gam + qnorm(0.975) * se.mod[1] )
  confgam <- confgam / pi * 180  
  confdel <- c(del - qnorm(0.975) * se.mod[2], del + qnorm(0.975) * se.mod[2] )
  confdel <- confdel / pi * 180
  confa <- c(a - qnorm(0.975) * se.a, a + qnorm(0.975) * se.a )
  confa <- confa / pi * 180
  confb <- c(b - qnorm(0.975) * se.b, b + qnorm(0.975) * se.b )
  confb <- confb / pi * 180
  confk <- c(k - qnorm(0.975) * se.k, k + qnorm(0.975) * se.k )
  info <- rbind( c(gam / pi * 180, confgam), c(del / pi * 180, confdel), 
  c(a / pi * 180, confa), c(b / pi * 180, confb),  c(k, confk) )
  rownames(info) <- c("gamma", "delta", "alpha", "beta", "kappa")
  colnames(info) <- c("estimate", "2.5%", "97.5%")
  modes <- rbind( c(a, b/2), c(a, b/2 + pi) )
  modes <- modes / pi * 180
  rownames(modes) <- c("mode 1", "mode 2")
  colnames(modes) <- c("co-latitude", "longitude")
  unitvectors <- cbind(m1, m2, m3)
  colnames(unitvectors) <- c("mu 1", "mu 2", "mu 3")
  list(info = info, modes = modes, unitvectors = unitvectors, loglik = ka$objective )
}