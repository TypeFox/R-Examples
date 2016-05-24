################################
#### Regularised discriminant analysis
#### Tsagris Michail 7/2012
#### mtsagris@yahoo.gr
#### References: Hastie T., Tibshirani R. & Friedman J. (2008)
#### Elements of Statistical Learning (2nd edition) p. 112-113. Springer
################################

rda <- function(xnew, x, ina, gam = 1, del = 0) {
  ## xnew is the new observation
  ## x contains the data
  ## ina is the grouping variable
  ## gam is between pooled covariance and diagonal
  ## gam*Spooled+(1-gam)*diagonal
  ## del is between QDA and LDA
  ## del*QDa+(1-del)*LDA
  ## mesi is NULL by default or it can be a matrix with the group means
  ## info can be NULL or it can be a list containing the
  ## covariance matrix of each group, the pooled covariance matrix
  ## and the spherical covariance matrix (this order must be followed)
  ## the mesi and info are particularly useful for the tuning of the rda, as
  ## they can speed the computations a lot.

  x <- as.matrix(x)
  n <- nrow(x)
  D <- ncol(x)
  xnew <- as.matrix(xnew)
  xnew <- matrix(xnew, ncol = D)
  nu <- nrow(xnew)  ## number of the new observations
  ina <- as.numeric(ina)
  nc <- max(ina)
  Ska <- array( dim = c(D, D, nc) )
  ta <- matrix(nrow = nu, ncol = nc)
  ng <- as.vector( table(ina) )
  ci <- log(ng / n)
  con <-  - D / 2 * log(2 * pi)  ## constant part
  sk <- array( dim = c(D, D, nc) )

  mesos <- aggregate( x, by = list(ina), mean )
  mesos <- as.matrix( mesos[, -1] )

  ni <- rep(ng - 1, each = D^2)

  for (m in 1:nc)  sk[, , m] <- cov( x[ina == m, ] )
  s <- ni * sk
  Sp <- apply(s, 1:2, sum) / (n - nc)  ## pooled covariance matrix
  sp <- diag( mean( diag(Sp) ), D ) ## spherical covariance matrix
  Sa <- gam * Sp + (1 - gam) * sp  ## regularised covariance matrix

  for (j in 1:nc) {
    Ska[, , j] <- del * sk[, , j] + (1 - del) * Sa
    ta[, j] <- ci[j] - 0.5 * log( det( Ska[, , j] ) ) -
      0.5 * mahalanobis( xnew, mesos[j, ], Ska[, , j] )
  }

  ta <- ta + con
  est <- apply(ta, 1, which.max)
  prob <- exp(ta) / rowSums( exp(ta) ) ## the probability of classification
  list(prob = prob, scores = ta, est = est)
}
