################################
#### ANOVA for hyper-spherical data (Log-likelihood ratio test)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 224
################################

lr.aov <- function(x, ina) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  x <- as.matrix(x)
  x <- x/sqrt(rowSums(x^2))  ## makes sure x are unit vectors
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  S <- aggregate(x, by = list(ina), sum)
  Ri <- sqrt( rowSums(S[, -1]^2) )  ## the resultant length of each group
  S <- colSums(x)
  R <- sqrt(sum(S^2))  ## the resultant length based on all the data
  ## Next we stimate the common concentration parameter kappa under H0 and H1
  Apk <- function(p, k) {
  besselI(k, p/2, expon.scaled = TRUE)/besselI(k, p/2 - 1, expon.scaled = TRUE)
  }
  Rk <- R/n
  k <- numeric(4)
  j <- 1
  k[j] <- Rk * (p - Rk^2)/(1 - Rk^2)
  j <- 2
  k[j] <- k[j - 1] - (Apk(p, k[j - 1]) - Rk)/(1 - Apk(p, k[j - 1])^2 -
  (p - 1)/k[j - 1] * Apk(p, k[j - 1]))
  while (abs(k[j] - k[j - 1]) > 1e-07) {
    j <- j + 1
    k[j] <- k[j - 1] - (Apk(p, k[j - 1]) - Rk)/(1 - Apk(p, k[j - 1])^2 -
    (p - 1)/k[j - 1] * Apk(p, k[j - 1]))
  }
  k0 <- k[j]  ## concentration parameter under H0
  Rk <- sum(Ri)/n
  k <- numeric(4)
  j <- 1
  k[j] <- Rk * (p - Rk^2)/(1 - Rk^2)
  j <- 2
  k[j] <- k[j - 1] - (Apk(p, k[j - 1]) - Rk)/(1 - Apk(p, k[j - 1])^2 -
  (p - 1)/k[j - 1] * Apk(p, k[j - 1]))
  while (abs(k[j] - k[j - 1]) > 1e-07) {
    j <- j + 1
    k[j] <- k[j - 1] - (Apk(p, k[j - 1]) - Rk)/(1 - Apk(p, k[j - 1])^2 -
    (p - 1)/k[j - 1] * Apk(p, k[j - 1]))
  }
  k1 <- k[j]  ## concentration parameter under H1
  apk0 <- (1 - p/2) * log(k0/2) + lgamma(p/2) +
  log(besselI(k0, p/2 - 1, expon.scaled = TRUE)) + k0
  apk1 <- (1 - p/2) * log(k1/2) + lgamma(p/2) +
  log(besselI(k1, p/2 - 1, expon.scaled = TRUE)) + k1
  w <- 2 * (k1 * sum(Ri) - k0 * R - n * apk1 + n * apk0)
  pvalue <- 1 - pchisq(w, (g - 1) * (p - 1))
  res <- c(w, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
