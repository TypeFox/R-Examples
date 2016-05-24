################################
#### Kuiper test of uniformity with circular data
#### Tsagris Michail 3/2016 
#### mtsagris@yahoo.gr
#### References: Jammalamadaka, S. Rao and SenGupta, A. (2001). 
#### Topics in Circular Statistics, pg. 153-155
################################

kuiper <- function(u, rads = FALSE, R = 1) {
  ## u is a vector with circular data
  ## if data are in rads set it to TRUE
  ## R is for Monte Carlo estimate of the p-value
  
  if (rads == FALSE)  u <- u / 180 * pi
  u <- sort(u) / (2 * pi)
  n <- length(u)
  i <- 1:n
  f <- sqrt(n)
  Vn <- f * ( max(u - (i - 1)/n) + max(i/n - u) )

  if (R == 1) {  ## asymptotic p-value is returned
    m2 <- (1:50)^2
    a1 <- 4 * m2 * Vn^2
    a2 <- exp(-2 * m2 * Vn^2)
    b1 <- 2 * ( a1 - 1 ) * a2
    b2 <- 8 * Vn / ( 3 * f ) * m2 * (a1 - 3) * a2
    pvalue <- sum(b1 - b2)

  } else {
    bvn <- numeric(R)
    for (j in 1:R) {
      x <- runif(n, 0, 2 * pi)
      x <- sort(x) / (2 * pi)
      bvn[j] <- f * ( max(x - (i - 1)/n) + max(i/n - x) )
    } 
    pvalue <- ( sum(bvn > Vn) + 1 ) / (R + 1)
  }

  res <- c(Vn, pvalue)
  names(res) <- c("Test", "p-value")
  res
}





