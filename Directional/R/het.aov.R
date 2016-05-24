################################
#### ANOVA for hyper-spherical data (Heterogeneous case, kappas not equal)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 228
################################

het.aov <- function(x, ina) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  x <- as.matrix(x)
  x <- x/sqrt(rowSums(x^2))  ## makes sure x are unit vectors
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  ni <- as.vector(table(ina))  ## group sample sizes
  kappa <- numeric(g)
    mi <- aggregate(x, by = list(ina), mean)
    mi <- mi[, -1]
    for (i in 1:g) {
    kappa[i] <- vmf(x[ina == i, ])$kappa
  }
  tw <- colSums(kappa * ni * mi)
  Tt <- 2 * ( sum(kappa * ni * sqrt(rowSums(mi^2))) - sqrt(sum(tw^2)) )
  pvalue <- 1 - pchisq(Tt, (p - 1) * (g - 1))
  res <- c(Tt, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
