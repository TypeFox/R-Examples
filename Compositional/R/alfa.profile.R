################################
#### Profile log-likelihood for choosing the value of alpha
#### Tsagris Michail 5/2013
#### References: Tsagris, M. T., Preston, S., and Wood, A. T. A. (2011).
#### A data-based power transformation for
#### compositional data. In Proceedings of the 4rth Compositional Data Analysis Workshop, Girona, Spain.
#### mtsagris@yahoo.gr
################################

alfa.profile <- function(x, a = seq(-1, 1, by = 0.01) ) {
  ## x contains the data
  ## a is the grid of values of the power parameter

  x <- as.matrix(x)  ## makes the data in a matrix form
  x <- x/rowSums(x)  ## makes sure the data are compositional
  D <- ncol(x)  ## number of components
  d <- D - 1  ## dimensionality of the simplex
  n <- nrow(x)  ## sample size of the data
  f <- (n - 1) / n
  qa <- numeric( length(a) )  ## the log-likelihood values will be stored here
  ja <- sum( log(x) )  ## part of the Jacobian of the alpha transformation
  con <-  - n/2 * d * log(2 * pi * f) - (n - 1) * d/2 + n * (d + 1/2) * log(D)

  for ( i in 1:length(a) ) {
    trans <- alfa( x, a[i] )
    aff <- trans$aff  ## the alpha-transformation
    sa <- trans$sa  ## part of the Jacobian determinant as well
    qa[i] <-  - n/2 * log( abs( det( cov(aff) ) ) ) + (a[i] - 1) * ja - D * sa
  }

  qa <- qa + con
  ## the green lines show a 95% CI for the true value of
  ## alpha using a chi-square distribution
  b <- max(qa) - qchisq(0.95, 1)/2
  plot(a, qa, type = "l", xlab = expression( paste(alpha, " values", sep = "") ),
  ylab = "Profile log-likelihood")
  abline(h = b, col = 2)
  ci <- c( min(a[qa >= b]), max(a[qa >= b]) )
  names(ci) <- paste(c("2.5", "97.5"), "%", sep = "")
  abline(v = ci[1], col = 3, lty = 2)
  abline(v = ci[2], col = 3, lty = 2)

  res <- c(a[which.max(qa)], max(qa), qa[a == 0])
  names(res) <- c('alfa', 'max.log.lik', 'log.lik0')
  list(result = res, ci = ci)
}
