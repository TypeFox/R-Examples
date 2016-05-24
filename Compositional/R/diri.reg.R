################################
#### Dirichlet regression for compositional data
#### Tsagris Michail 6/2014
#### mtsagris@yahoo.gr
################################

diri.reg <- function(y, x, plot = TRUE, xnew = NULL) {
  ## y is the compositional data
  y <- as.matrix(y)
  n <- nrow(y)
  y <- y/rowSums(y)
  x <- as.matrix( cbind(1, x) )
  ## the line above makes sure y is compositional data and
  ## then the unit vector is added to the desing matrix
  d <- ncol(y) - 1  ## dimensionality of the simplex
  z <- list(y = log(y), x = x)

    dirireg <- function(param, z = z) {
      ## param contains the parameter values
      ## z contains the compositional data and independent variable(s)
      phi <- exp( param[1] )  ## this avoids negative values in phi
      para <- param[-1]
      y <- z$y
      x <- z$x
      ## y is the compositional data and xa the independent variable(s)
      n <- nrow(y)  ## sample size
      d <- ncol(y) - 1  ## dimensionality of  the simplex
      be <- matrix(para, ncol = d)  ## puts the beta parameters in a matrix
      mu1 <- cbind(1, exp(x %*% be))
      ma <- mu1/rowSums(mu1)  ## the fitted values
      l <- -( n * lgamma(phi) - sum( lgamma(phi * ma) ) +
      sum( diag( y %*% t(phi * ma - 1) ) )  )
      ## l is the log-likelihood
    l
  }

  runtime <- proc.time()
  rla <- log(y[, -1] / y[, 1])  ## additive log-ratio transformation
  ini <- as.vector( coef(lm.fit(x, rla)) )  ## initial values
  ## based on the logistic normal
  ## the next lines optimize the dirireg function and
  ## estimate the parameter values

  el <- NULL
  options(warn = -1)
  qa <- nlm(dirireg, c(3, as.vector( t(ini)) ), z = z)
  el[1] <- -qa$minimum
  qa <- nlm(dirireg, qa$estimate, z = z)
  el[2] <- -qa$minimum
  vim <- 2
  while (el[vim] - el[vim - 1] > 1e-06) {
    ## the tolerance value can of course change
    vim <- vim + 1
    qa <- nlm(dirireg, qa$estimate, z = z)
    el[vim] <- -qa$minimum
  }
  qa <- nlm(dirireg, qa$estimate, z = z, hessian = TRUE)
  log.phi <- qa$estimate[1]
  para <- qa$estimate[-1]  ## estimated parameter values
  beta <- matrix(para, ncol = d)  ## matrix of the betas
  colnames(beta) <- colnames(y[, -1])  ## names of the betas
  seb <- sqrt( diag( solve(qa$hessian) ) )  ## std of the estimated betas
  std.logphi <- seb[1]  ## std of the estimated log of phi
  seb <- matrix(seb[-1], ncol = d)  ## std of the estimated betas

  if ( !is.null( colnames(y) ) ) {
    colnames(seb) <- colnames(y[, -1])
  } else  colnames(seb) <- paste("Y", 1:d, sep = "")

  if ( !is.null(xnew) ) {
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    mu <- cbind( 1, exp(xnew %*% beta) )
    est <- mu / rowSums(mu)
  } else {
    mu <- cbind( 1, exp(x %*% beta) )
    est <- mu / rowSums(mu)  ## fitted values
    lev <- ( exp(log.phi) + 1 ) * rowSums( (y - est)^2 / mu )
    if (plot == TRUE) {
      plot(1:n, lev, main = "Influence values", xlab = "Observations",
      ylab = expression( paste("Pearson ", chi^2, "statistic") ) )
      lines(1:n, lev, type = "h")
      abline(h = qchisq(0.95, d), lty = 2, col = 2)
    }
  }

  runtime <- proc.time() - runtime

  if ( is.null(colnames(x)) ) {
    p <- ncol(x) - 1
    rownames(beta) <- c("constant", paste("X", 1:p, sep = "") )
    if ( !is.null(seb) )  rownames(seb) <- c("constant", paste("X", 1:p, sep = "") )
  } else {
    rownames(beta)  <- c("constant", colnames(x)[-1] )
    if  ( !is.null(seb) ) rownames(seb) <- c("constant", colnames(x)[-1] )
  }

  list(runtime = runtime, loglik = -qa$minimum, phi = exp(log.phi), log.phi = log.phi,
  std.logphi = std.logphi, beta = beta, seb = seb, lev = lev, est = est)
}
