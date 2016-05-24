################################
#### Estimating parameters for compositional data
#### using the additive or the isometric log-ratio transformation
#### Tsagris Michail 1/2016
#### mtsagris@yahoo.gr
################################

comp.den <- function(x, type = "alr", dist = "normal") {
  ## x is the compositional data
  ## type is the type of transformation, "alr", "ilr"
  ## dist is the distribution to be fitted,
  ## "normal", "rob", "spatial", "t", "skewnorm"
  x <- as.matrix(x)  ## makes sure x is a matrix
  x <- x / rowSums(x)  ## makes sure x is compositional data
  ## if type = "none" or dist = "dirichlet" the Dirichlet is fitted
  if (dist == "normal") {
    if (type == "alr") {  ## additive log-ratio transformation
      y <- log(x[, -1]/ x[, 1])
      m <- colMeans(y)
      mu <- c( 1, exp(m) )
      mu <- mu / sum(mu)
      s <- cov(y)
    } else {
      y <- alfa(x, 0)
      m <- colMeans(y)
      mu <- alfainv(m, 0)
      s <- cov(y)
    }
    result <- list(mean = m, comp.mean = mu, covariance = s)
  } else if (dist == "t") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- log(x[, -1]/ x[, 1])
        mod <- multivt(y)
        m <- mod$center
        mu <- c( 1, exp(m) )
        mu <- mu / sum(mu)
        s <- mod$scatter
        dof <- mod$dof
      } else {
        y <- alfa(x, 0)
        mod <- multivt(y)
        m <- mod$center
        mu <- alfainv(m, 0)
        s <- mod$scatter
        dof <- mod$dof
      }
      result <- list(mean = m, comp.mean = mu, covariance = s, dof = dof)
  } else if (dist == "rob") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- log(x[, -1]/ x[, 1])
        mod <- MASS::cov.rob(y, method = "mcd")
        m <- mod$center
        mu <- c( 1, exp(m) )
        mu <- mu / sum(mu)
        s <- mod$cov
        best <- mod$best
      } else {
        y <- alfa(x, 0)
        mod <- MASS::cov.rob(y, method = "mcd")
        m <- mod$center
        mu <- alfainv(m, 0)
        s <- mod$cov
        best <- mod$best
      }
      result <- list(mean = m, comp.mean = mu, covariance = s, best = best)
  } else if (dist == "spatial") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- log(x[, -1]/ x[, 1])
        delta <- spat.med(y)
        comp.delta <- c( 1, exp( delta ) )
        comp.delta <- delta / sum( delta )
        s <- sscov(y)
      } else {
        y <- alfa(x, 0)
        delta <- spat.med(y)
        comp.delta <- alfainv(delta, 0)
        s <- sscov(y)
      }
      result <- list(spat.med = delta, comp.spat.med = comp.delta, ssc = s)
  } else if (dist == "skewnorm") {
      if (type == "alr") {  ## additive log-ratio transformation
        y <- log(x[, -1]/ x[, 1])
        mod <- sn:: msn.mle(y = y)
        beta <- as.vector( mod$dp$beta )
        Omega <- as.matrix( mod$dp$Omega )
        alpha <- as.vector(mod$dp$alpha)
        cobeta <- c( 1, exp( beta) )
        cobeta <- beta / sum(beta)
      } else {
        y <- alfa(x, 0)
        mod <- sn:: msn.mle(y = y)
        beta <- as.vector( mod$dp$beta )
        Omega <- as.matrix( mod$dp$Omega )
        alpha <- as.vector(mod$dp$alpha)
        cobeta <- alfainv(beta, 0)
      }
      result <- list(beta = beta, Omega = Omega, alpha = alpha, comp.beta = cobeta)
  }
  result
}
