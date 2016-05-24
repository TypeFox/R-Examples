L_hNV <- function(p, y = y, X = X, sc = sc) {
  
  b <- p[3:length(p)]
  sigmau2_dmu <- p[1]
  sigmav2 <- p[2]
  
  N <- length(y)
  
  epsilon <- (y - X%*%b)
  
  ret = (N * log(sqrt(2) / sqrt(pi)) + N * log(1 / (sqrt(sigmau2_dmu + sigmav2)))
         + sum(log(pnorm(-sc*(epsilon * (sqrt(sigmau2_dmu) / sqrt(sigmav2))) / 
                           (sqrt(sigmau2_dmu + sigmav2)))))
         - 1 / (2 * (sigmau2_dmu + sigmav2)) * sum(epsilon^2))
  
  names(ret) <- "Log-Lik normal/half-normal distribution (rho=0)"
  return(ret)
  
}
