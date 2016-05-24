################################
#### ANOVA for hyper-spherical data (Embedding approach)
#### Tsagris Michail 1/2015
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics, page 225-226
################################

embed.aov <- function(x, ina) {
  ## x contains all the data
  ## ina is an indicator variable of each sample
  ina <- as.numeric(ina)
  g <- max(ina)  ## how many groups are there
  ni <- as.vector(table(ina))
  x <- as.matrix(x)
  x <- x/sqrt( rowSums(x^2) )  ## makes sure x are unit vectors
  p <- ncol(x)  ## dimensionality of the data
  n <- nrow(x)  ## sample size of the data
  S <- aggregate(x, by = list(ina), mean)
  Rbi <- sqrt( rowSums(S[, -1]^2) ) ## the mean resultant length of each group
  S <- colMeans(x)
  Rbar <- sqrt(sum(S^2))  ## the mean resultant length based on all the data
  Ft <- ((n - g) * (p - 1) * (sum(ni * Rbi^2) - n * Rbar^2))/((g - 1) * (p - 1) *
  (n - sum(ni * Rbi^2)))
  pvalue <- 1 - pf(Ft, (g - 1) * (p - 1), (n - g) * (p - 1))
  res <- c(Ft, pvalue)
  names(res) <- c('test', 'p-value')
  res
}
