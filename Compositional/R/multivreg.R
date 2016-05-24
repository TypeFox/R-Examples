################################
#### Multivariate regression
#### Tsagris Michail 6/2011
#### mtsagris@yahoo.gr
#### References: Mardia K.V., kent J.T. & Bibby J.M. (1979)
#### Multivariate Analysis p. 318-320. Academic Press
################################

multivreg <- function(y, x, plot = TRUE, xnew = NULL) {
  ## y is the dependent variable and must be a matrix
  ## with at least two columns
  ## x contains the independent variable(s) which have to be
  ## in a matrix format or a vector if you have just one
  y <- as.matrix(y)
  x <- as.matrix(x)
  n <- nrow(y)  ## sample size
  d <- ncol(y)  ## dimensionality of y
  p <- ncol(x)  ## dimensionality of x
  mod <- lm(y ~ x)   ## linear regression
  res <- resid(mod)  ## residuals
  s <- cov(res) * (n - 1) / (n - p - 1)
  sxx <- cov(x)  ## covariance of the independent variables
  dres <- sqrt( mahalanobis(res, numeric(d), s) )  ## Mahalanobis distances
  ## of the residuals
  mx <- colMeans(x)  ## mean vector of the independent variales
  dx <- sqrt( mahalanobis(x, mx, sxx) )  ## Mahalanobis distances
  ## of the independent variables
  crit.res <- sqrt( qchisq(0.975, d) )
  crit.x <- sqrt( qchisq(0.975, p) )

  if (plot == TRUE) {
    plot(dx, dres, xlim = c(0, max(dx) + 0.5), ylim = c(0, max(dres) + 0.5),
    xlab = "Mahalanobis distance of x", ylab = "Mahalanobis distance
    of residuals")
    abline(h = crit.res)
    abline(v = crit.x)
  }

  resid.out <- as.vector( which(dres > crit.res) )
  x.leverage <- which(dx > crit.x)
  out.and.lever <- which(dx > crit.x & dres > crit.res)

  if ( is.null(xnew) ) {
    est <- fitted(mod)
  } else {
    xnew <- cbind(1, xnew)
    xnew <- as.matrix(xnew)
    est <- xnew %*% coef(mod)
  }

  moda <- summary(mod)
  suma <- array(dim = c(1 + p, 6, d))
  r.squared <- numeric(d)
  mse <- deviance(mod)/( n - p - 1 )

  for (i in 1:d) {
    wa <- as.matrix( coef(moda[[i]]) )
    wa <- cbind( wa, wa[, 1] - qt(0.975, n - p - 1) * mse[i] ,
    wa[, 1] + qt(0.975, n - p - 1) * mse[i] )
    colnames(wa)[5:6] <- paste(c(2.5, 97.5), "%", sep = "")
    suma[, , i] <- wa
    r.squared[i] <- as.numeric( moda[[i]]$r.squared )
  }

  if ( is.null(colnames(y)) ) {
    dimnames(suma) <- list( rownames(wa), colnames(wa),
    paste("Y", 1:d, sep = "") )
    names(r.squared) <- paste("Y", 1:d, sep = "")
    colnames(est) <- paste("Y", 1:d, sep = "")
  } else {
    dimnames(suma) <- list( rownames(wa), colnames(wa), colnames(y) )
    names(r.squared) <- colnames(y)
    colnames(est) <- colnames(y)
  }

  list(suma = suma, r.squared = r.squared, resid.out  = resid.out,
  x.leverage = x.leverage, out = out.and.lever, est = est)
}

