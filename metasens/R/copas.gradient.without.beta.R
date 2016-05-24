copas.gradient.without.beta <- function(x, gamma, TE, seTE){
  
  
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
  ## this is towards the smalles value for log
  ##
  v[v < -37] <- -37
  ##
  ##
  ci2 <- lambda(u)*(u+lambda(u))
  
  
  ##
  ## Derivatives for mu:
  ##
  ## term 1 is zero
  ##
  ## term 2:
  ##
  grad.mu <- (TE-mu)/(tau^2+sigma^2)
  ##
  ## term 3 is always 0
  ##
  ## term 4:
  ##
  grad.mu <- (grad.mu -
              (rho.tilde/(sqrt((tau^2+sigma^2)*(1-rho.tilde^2))))*lambda(v)
              )
  
  ##
  ## Derivatives for rho:
  ##
  ## term 1:
  ##
  grad.rho <- (-ci2*rho*sigma^2*sigma^2/
               (seTE^2*(tau^2+sigma^2)))
  ##
  ## term 2:
  ##
  grad.rho <- (grad.rho +
               (((TE-mu)^2)*
                ci2*rho*sigma^2*sigma^2/
                (seTE^2*((tau^2+sigma^2)^2))
                )
               )
  ##
  ## term 4:
  ##
  top <- u + rho.tilde*(TE-mu)/sqrt(tau^2+sigma^2) 
  bottom <- sqrt((1-rho.tilde^2))
  ##
  diff.top <- ((top-u)/rho - (top-u)*rho*ci2/(1-ci2*rho^2) +
               2*(top-u)*rho*tau^2*ci2/(seTE^2+tau^2*(1-ci2*rho^2))
               )
  ##
  eta <- seTE^2/(seTE^2 + tau^2*(1-ci2*rho^2))
  ##
  diff.bottom <- ((1-rho.tilde^2)^(-1.5))*rho*eta*(1 + (rho^2*tau^2*ci2*eta)/seTE^2 )
  ##
  grad.rho <- grad.rho + (top*diff.bottom+diff.top/bottom)*lambda(v)
  
  ##
  ## gradient for square-root of variance (tau)
  ##
  ## term 1:
  ##
  grad.tau <- -0.5/(tau^2+sigma^2)
  ##
  ## term 2:
  ##
  grad.tau <- (grad.tau +
               0.5*((TE-mu)^2)/((tau^2+sigma^2)^2)
               )
  ##
  ## term 4:
  ##
  grad.tau <- (grad.tau +
               (-sqrt(sigma^2)*
                (TE-mu)*
                rho/(bottom*((tau^2+sigma^2)^2)) -
                0.5*top*((1-rho.tilde^2)^(-1.5))*
                sigma^2*rho^2/((tau^2+sigma^2)^2)
                )*lambda(v)
               )
  ##
  grad.tau <- 2 * tau * grad.tau
  
  
  ##  
  ## negative gradient as seek to minimise - log likelihood
  ##
  res <- c(sum(-grad.mu), sum(-grad.rho), sum(-grad.tau))
  ##
  res
}
