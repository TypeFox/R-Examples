################################
#### Permutation based hypothesis testing 
#### for a zero correlation coefficient 
####
################################

permcor <- function(x, R = 999) {
  ## x is a 2 column matrix containing the data
  ## type can be either "pearson" or "spearman"
  ## R is the number of permutations
  
  x <- as.matrix(x)
  n <- nrow(x)
  r <- cor(x)[2]
  test <- 0.5 * log( (1 + r)/(1 - r) )  ## the test statistic
  x2 <- x[, 2]
  m1 <- mean(x[, 1])    ;   m12 <- sum(x[, 1]^2)
  m2 <- mean(x2)        ;   m22 <- sum(x2^2)
  
  sxy <- numeric(R)
  for (i in 1:R) {
    x1 <- sample(x[ , 1], n)
    sxy[i] <- sum(x1 * x2)
  }  
    rb <- (sxy - n * m1 * m2) / sqrt( (m12 - n * m1^2) * (m22 - n * m2^2) )
    tb <- 0.5 * log( (1 + rb)/(1 - rb) )  ## the test statistic

  pvalue <- ( sum( abs(tb) > abs(test) ) + 1 ) / (R + 1)  ## bootstrap p-value
  res <- c( r, pvalue )
  names(res) <- c('correlation', 'p-value')
  res
}


