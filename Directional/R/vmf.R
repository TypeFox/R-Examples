################################
#### MLE for a von-Mises Fisher distribution 
#### Tsagris Michail 10/2013 
#### mtsagris@yahoo.gr
#### References: Suvrit Sra (2012, Computational Statistics)
#### A short note on parameter approximation for von Mises-Fisher distributions: 
#### and a fast implementation of Is(x)
################################

vmf <- function(x, tol = 1e-06) {
  ## x contains the data
  ## tol specifies the tolerance value for convergence
  ## when estimating the concentration parameter
  x <- as.matrix(x)
  x <- x/sqrt(rowSums(x^2))
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  Apk <- function(p, k) {
   besselI(k, p/2, expon.scaled = T)/besselI(k, p/2 - 1, expon.scaled = T)
  }
  m1 <- colSums(x)
  R <- sqrt(sum(m1^2))/n  ## mean resultant length
  m <- m1/(n * R)
  k <- numeric(4)
  i <- 1
  k[i] <- R * (p - R^2)/(1 - R^2)
  if (k[i] > 100000) { 
    k <- k[i] 
  } else { 
    i <- 2
    k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 - 
    (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
    while (abs(k[i] - k[i - 1]) > tol) {
      i <- i + 1
      k[i] <- k[i - 1] - (Apk(p, k[i - 1]) - R)/(1 - Apk(p, k[i - 1])^2 - 
      (p - 1)/k[i - 1] * Apk(p, k[i - 1]))
    }  
    k <- k[i]
  }
  loglik <- n * (p/2 - 1) * log(k) - 0.5 * n * p * log(2 * pi) - 
  n * (log(besselI(k, p/2 - 1, expon.scaled = T)) + k) + k * sum(x %*% m)
  vark <- 1/(n * (1 - Apk(p, k)/k - Apk(p, k)^2))
  list(mu = m, kappa = k, MRL = R, vark = vark, loglik = loglik)
}