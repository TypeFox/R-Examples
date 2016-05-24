loglik.2binormal <-
function(z.1, z.2, mu, sigma, rho, p){

  l.m <- sum(d.binormal(z.1, z.2, 0, 1, 0)+log(p*exp(d.binormal(z.1, z.2, mu, sigma, rho)-d.binormal(z.1, z.2, 0, 1, 0))+(1-p)))
  
  return(l.m) 
}

