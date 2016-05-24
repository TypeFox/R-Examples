################################
#### Discrminant analysis for directional data
#### assuming a von Mises-Fisher distribution
#### Cross-validation for the performance
#### Tsagris Michail 03/2014
#### mtsagris@yahoo.gr
#### References: J. E. Morris and P. J. Laycock (1974)
#### Discriminant Analysis of Directional Data (Biometrika)
################################

vmf.da <- function(x, ina, fraction = 0.2, R = 1000, seed = FALSE) {
  ## x is the data set
  ## ina is the group indicator variable
  ## fraction denotes the percentage of the sample to be used as the test sample
  ## R is the number of cross validations
  x <- as.matrix(x)
  x <- x / sqrt(rowSums(x^2))  ## makes sure x is unit vectors
  p <- ncol(x)  ## p is the dimensionality of the data
  per <- numeric(R)
  n <- nrow(x)  ## sample size
  ina <- as.numeric(ina)
  frac <- round(fraction * n)
  g <- max(ina)  ## how many groups are there
  mesi <- matrix(nrow = g, ncol = p)
  k <- numeric(g)
  ## if seed==TRUE then the results will always be the same
  if (seed == TRUE)  set.seed(1234567)
  for (i in 1:R) {
    mat <- matrix(nrow = frac, ncol = g)
    est <- numeric(frac)
    nu <- sample(1:n, frac)
    test <- x[nu, ]
    id <- ina[-nu]
    train <- x[-nu, ]
    for (j in 1:g) {
      da <- vmf(train[id == j, ])  ## estimates the parameters of the vMF
      mesi[j, ] <- da$mu  ## mean direction of each group
      k[j] <- da$kappa ## concentration of each group
    }
    for (j in 1:g) {
      mat[, j] <- (p/2 - 1) * log(k[j]) + k[j] * test %*% mesi[j, ] - 0.5 *
        p * log(2 * pi) - log(besselI(k[j], p/2 - 1, expon.scaled = TRUE)) - k[j]
    }
    est <- apply(mat, 1, which.max)
    per[i] <- sum(est == ina[nu])/frac
  }
  percent <- mean(per)
  s1 <- sd(per)
  s2 <- sqrt(percent * (1 - percent)/R)
  conf1 <- c(percent - 1.96 * s1, percent + 1.96 * s1)  ## 1st way of a CI
  conf2 <- c(percent - 1.96 * s2, percent + 1.96 * s2)  ## 2nd way of a CI
  ## next we check if the confidence limits exceeds the allowed limits
  if (conf1[2] > 1) conf1[2] <- 1
  if (conf1[1] < 0) conf1[1] <- 0
  if (conf2[2] > 1) conf2[2] <- 1
  if (conf2[1] < 0) conf2[1] <- 0
  conf3 <- quantile(per, probs = c(0.025, 0.975))  ## 3rd way of a CI
  ci <- rbind(conf1, conf2, conf3)
  colnames(ci) <- c("2.5%", "97.5%")
  rownames(ci) <- c("standard", "binomial", "empirical")
  percent <- c(percent, s1, s2)
  names(percent) <- c('percent', 'sd1', 'sd2')
  list(percent = percent, ci = ci)
}
