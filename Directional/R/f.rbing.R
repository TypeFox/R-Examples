################################
#### Simulating from a Bingham distribution
#### Tsagris Michail 02/2014
#### mtsagris@yahoo.gr
#### References: A new method to simulate the Bingham and related distributions
#### in directional data analysis with applications
#### Kent J.T., Ganeiber A.M. and Mardia K.V. (2013)
#### http://arxiv.org/pdf/1310.8110v1.pdf
#### and
#### Exact Bayesian Inference for the Bingham Distribution
#### C.J. Fallaize and T. Kypraios (2014)
#### http://arxiv.org/pdf/1401.2894v1.pdf
################################

######### Christopher Fallaize and Theodore Kypraios' code

f.rbing <- function(n, lam) {
  ## n is the sample size
  ## lam are the q - 1 non zero eigenvalues

  lam <- sort(lam, decreasing = TRUE)  ## sort the eigenvalues in desceding order
  nsamp <- 0
  X <- NULL
  lam.full <- c(lam, 0)
  qa <- length(lam.full)
  mu <- numeric(qa)
  A <- diag(lam.full)
  SigACG.inv <- diag(qa) + 2 * A
  SigACG <- solve(SigACG.inv)
  Ntry <- 0

  while (nsamp < n) {
    x.samp <- FALSE
    while (x.samp == FALSE) {
      yp <- MASS::mvrnorm(n = 1, mu = mu, Sigma = SigACG)
      y <- yp / sqrt( t(yp) %*% yp )
      lratio <-  - mahalanobis(y, mu, A, inverted = TRUE) - qa/2 * log(qa) +
      0.5 * (qa - 1) + qa/2 * log( mahalanobis( y, mu, SigACG.inv, inverted = TRUE ) )
      if ( log(runif(1) ) < lratio) {
        X <- c(X, y)
        x.samp <- TRUE
        nsamp <- nsamp + 1
      }
      Ntry <- Ntry + 1
    }
  }

  if (n > 1)
    X <- matrix(X, byrow = TRUE, ncol = qa)
  ## the X contains the simulated values
  ## the avtry is the estimate of the M in rejection sampling
  ## 1/M is the probability of acceptance
  list(X = X, avtry = Ntry/n)
}
