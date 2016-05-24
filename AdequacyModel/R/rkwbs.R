rkwbs <- function(n, alpha, beta, a, b){
  
  V = (1 - (1-runif(n,0,1))^(1/b))^(1/a)
  InversePhi = qnorm(V,0,1)
  return(beta*((alpha*InversePhi)/2 + (1+(alpha^(2))*(InversePhi^2))^(0.5))^(2) )
}