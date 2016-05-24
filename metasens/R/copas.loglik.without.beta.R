copas.loglik.without.beta <- function(x, gamma=c(-1.5,0.08),
                                      TE, seTE){
  
  
  mu  <- x[1]
  rho <- x[2]
  tau <- x[3]
  ##
  ## TE   <=> estimated treatment effect
  ## seTE <=> standard error from trials, conditional on publication
  
  
  ## Copas, Shi (2000), Biostatistics, p. 250:
  ##
  u <- gamma[1] + gamma[2]/seTE
  ##
  sigma <- sqrt(seTE^2/(1-rho^2*lambda(u)*(u+lambda(u))))
  rho.tilde <- rho*sigma/sqrt(tau^2+sigma^2)
  ##
  v <- ((u +
         rho.tilde *
         (TE-mu)/(sqrt(tau^2+sigma^2))
         ) /
        sqrt(1-rho.tilde^2)
        )
  ##
  ## avoid numerical problems by replacing 0's in pnorm(v):
  ## qnorm(1e-320) = -38.26913
  ## this is towards the smallest value for log
  ##
  v[v < -37] <- -37
  ##
  ## take minus log-likelihood and minimise it;
  ## leave out log(pnorm(u)) as this is a constant
  ##
  ell <- -(-0.5*log(tau^2+sigma^2) -
           (TE-mu)^2 / (2*(tau^2+sigma^2)) +
           log(pnorm(v))
           )
  
  res <- sum(ell)
  ##
  res
}
