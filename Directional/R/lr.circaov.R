################################
#### ANOVA for cicular data (Likelihood ratio test)
#### Tsagris Michail 1/2015
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 136
################################

lr.circaov <- function(u, ina, rads = FALSE) {
  ## u contains all the circular data in radians or degrees
  ## ina is an indicator variable of each sample
  n <- length(u)  ## sample size
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ## if the data are in degrees we transform them into radians
  if (rads == F)  u <- u * pi/180
  x <- cbind(cos(u), sin(u))
  Ci <- aggregate(x[, 1], by = list(ina), sum)
  Si <- aggregate(x[, 2], by = list(ina), sum)
  Ri <- sqrt(Ci[, 2]^2 + Si[, 2]^2)  ## the resultant length of each group
  mi <- aggregate(x, by = list(ina), mean)
  mi <- mi[, -1]
  mi <- mi/sqrt(rowSums(mi^2))  ## mean direction of each group
  m <- colMeans(x)
  m <- m/sqrt(sum(m^2))  ## mean direction based on all the data
  m <- matrix(rep(m, g), nrow = g, byrow = TRUE)
  C <- sum( cos(u) )
  S <- sum( sin(u) )
  R <- sqrt(C^2 + S^2)  ## the resultant length based on all the data
  ## Next we estimate the common concentration parameter kappa
  kappa <- circ.summary(u, rads = TRUE, plot = FALSE)$kappa
  ## kappa is the estimated concentration parameter based on all the data
  w <- kappa * sum(Ri * rowSums((mi - m)^2))
  pvalue <- 1 - pchisq(w, g - 1)
  res <- c(w, pvalue, kappa)
  names(res) <- c('test', 'p-value', 'kappa')
  res
}
