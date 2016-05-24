###############################
#### Linear cicrular correlation
#### Tsagris Michail 3/2014
#### mtsagris@yahoo.gr
#### References: Mardia Kanti V. and Jupp Peter E. (2000)
#### Directional statistics
################################

circlin.cor <- function(theta, x, rads = FALSE) {
  ## theta is a angular variable in degrees by default
  ## x is euclidean variable or a matrix containing euclidean variables
  x <- as.matrix(x)
  n <- length(theta)  ## sample size
  if (rads == FALSE)  theta <- theta * pi/180
  rxc <- as.numeric( cor(cos(theta), x) )  ## and cos(theta) and x correlation
  rxs <- as.numeric( cor(sin(theta), x) )  ## sin(theta) and x correlation
  rcs <- cor(cos(theta), sin(theta))  ## cos(theta) and sin(theta) correlation
  R2xt <- (rxc^2 + rxs^2 - 2 * rxc * rxs * rcs)/(1 - rcs^2)
  ## linear-circular correlation
  Ft <- (n - 3) * R2xt/(1 - R2xt)  ## F-test statistic value
  pvalue <- 1 - pf(Ft, 2, n - 3)
  res <- cbind(R2xt, pvalue)
  colnames(res) <- c('R-squared', 'p-value')
  res
}
