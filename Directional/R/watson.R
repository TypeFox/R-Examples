################################
#### Watson test of uniformity with circular data
#### Tsagris Michail 3/2016
#### mtsagris@yahoo.gr
#### References: Jammalamadaka, S. Rao and SenGupta, A. (2001).
#### Topics in Circular Statistics, pg. 156-157
################################

watson <- function(u, rads = FALSE, R = 1) {
  ## u is a vector with circular data
  ## if data are in rads set it to TRUE
  ## R is for Monte Carlo estimate of the p-value

  if (rads == FALSE)  u <- u / 180 * pi
  u <- sort(u) / (2 * pi)
  n <- length(u)
  i <- 1:n
  Wn <- sum( ( ( u - (i - 0.5)/n ) - ( mean(u) - 0.5 ) )^2 ) + 1 / ( 12 * n )

  if (R == 1) {  ## asymptotic p-value is returned
    m <- 1:20
    pvalue <- 2 * sum( ( - 1 )^( m - 1 ) * exp(-2 * m^2 * pi^2 * Wn) )

  } else {
    bwn <- numeric(R)
    for (j in 1:R) {
      x <- runif(n, 0, 2 * pi)
      x <- sort(x) / (2 * pi)
      bwn[j] <- ( max(x - (i - 1)/n) + max(i/n - x) )
    }
    pvalue <- ( sum(bwn > Wn) + 1 ) / (R + 1)
  }

  res <- c(Wn, pvalue)
  names(res) <- c("Test", "p-value")
  res
}

