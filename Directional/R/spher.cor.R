################################
#### Spherical-spherical correlation 
#### Tsagris Michail 4/2014 
#### mtsagris@yahoo.gr
#### References: Kanti V. Mardia and Peter E. Jupp
#### Directional statistics p.g. 254-255 
################################

spher.cor <- function(x, y) {
  ## x and y are two (hyper-)spherical variables
  x <- as.matrix(x)
  y <- as.matrix(y)
  x <- x/sqrt(rowSums(x^2))
  y <- y/sqrt(rowSums(y^2))
  stand <- function(x) x - mean(x)
  p <- ncol(x)  ## dimension of x
  q <- ncol(y)  ## dimension of y
  x <- apply(x, 2, stand)  ## subtract the mean
  y <- apply(y, 2, stand)  ## subtract the mean
  n <- nrow(x)  ## sample size
  s11 <- crossprod(x) / n
  s12 <- crossprod( x, y ) / n
  s21 <- t( s12 )
  s22 <- crossprod(y) / n
  rsq <- sum( diag(solve(s11, s12) %*% solve(s22, s21)) )
  test <- n * rsq
  pvalue <- 1 - pchisq(test, p * q)
  res <- c(rsq, pvalue)
  names(res) <- c('R-squared', 'p-value')
  res
}