# This function is modified from HWE.exact {genetics}

HWE.exact <- function (g){
  
  if(!all(g %in% 0:2)){
    return(1)
  }
  
  n11 <- sum(g == 0, na.rm = TRUE)
  n12 <- sum(g == 1, na.rm = TRUE)
  n22 <- sum(g == 2, na.rm = TRUE)
  
  n1 <- 2 * n11 + n12
  n2 <- 2 * n22 + n12
  dhwe2 <- function(n11, n12, n22) {
    f <- function(x) lgamma(x + 1)
    n <- n11 + n12 + n22
    n1 <- 2 * n11 + n12
    n2 <- 2 * n22 + n12
    exp(log(2) * (n12) + f(n) - f(n11) - f(n12) - f(n22) - 
          f(2 * n) + f(n1) + f(n2))
  }
  x12 <- seq(n1%%2, min(n1, n2), 2)
  x11 <- (n1 - x12)/2
  x22 <- (n2 - x12)/2
  dist <- data.frame(n11 = x11, n12 = x12, n22 = x22, density = dhwe2(x11, 
                                                                      x12, x22))
  dist <- dist[order(dist$density), ]
  STATISTIC <- c(N11 = n11, N12 = n12, N22 = n22)
  PARAMETER <- c(N1 = n1, N2 = n2)
  PVAL <- cumsum(dist$density)[dist$n11 == n11 & dist$n12 == 
                                 n12 & dist$n22 == n22]
  PVAL
}
