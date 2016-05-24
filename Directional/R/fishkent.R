################################
#### H0: Fisher versus H1: Kent distribution
#### Tsagris Michail 06/2014 
#### mtsagris@yahoo.gr
#### References: Louis-Paul Rivest (1986)
#### Statistics & Probability Letters, 4: 1-4.
#### Modified Kent's statistics for testing goodness of fit for the 
#### Fisher distribution in small concentrated samples
################################

fishkent <- function(x, B = 999) {
  ## x contains the data
  ## B is by default eaual to 999 bootstrap re-samples
  ## If B==1 then no bootstrap is performed
  n <- nrow(x)  ## sample size
  estim <- vmf(x)
  k <- estim$k  ## the estimated concentration parameter
  ## under the H0, that the Fisher distribution is true
  mu <- estim$mu  ## the estimated mean direction under H0
  e1 <- c(1, 0, 0)
  i3 <- diag(3)
  P <- i3 -  tcrossprod(e1 - mu) / (1 - mu[1])
  y <- tcrossprod(x, P)[, 2:3]
  lam <- eigen( crossprod(y) / n )$values
  rat <- besselI(k, 0.5, expon.scaled = T)/besselI(k, 2.5, expon.scaled = T)
  T <- n * (k/2)^2 * rat * (lam[1] - lam[2])^2
  if (B == 1) {
    pvalue <- 1 - pchisq(T, 2)
    res <- c(T, pvalue)
    names(res) <- c('test', 'p-value')
  } else {
    Tb <- numeric(B)
    for (i in 1:B) {
      nu <- sample(1:n, n, replace = T)
      estim <- vmf(x[nu, ])
      k <- estim$k  ## the estimated concentration parameter
      ## under the H0, that the Fisher distribution is true
      mu <- estim$mu  ## the estimated mean direction under H0
      P <- i3 -  tcrossprod(e1 - mu) / (1 - mu[1])
      y <- tcrossprod(x[nu, ], P)[, 2:3]
      lam <- eigen( crossprod(y) / n )$values
      rat <- besselI(k, 0.5, expon.scaled = T)/besselI(k, 2.5, expon.scaled = T)
      Tb[i] <- n * (k/2)^2 * rat * (lam[1] - lam[2])^2
    }
    res <- c( T, (sum(Tb > T) + 1)/(B + 1) )
    names(res) <- c('test', 'Bootstrap p-value')
  }
  res 
}