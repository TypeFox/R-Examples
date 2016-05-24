d.binormal <-
function(z.1, z.2, mu, sigma, rho){

  loglik <- (-log(2)-log(pi)-2*log(sigma) - log(1-rho^2)/2 - (0.5/(1-rho^2)/sigma^2)*((z.1-mu)^2 -2*rho*(z.1-mu)*(z.2-mu) + (z.2-mu)^2))

  return(loglik)
}

