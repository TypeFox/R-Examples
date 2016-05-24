size_n.regII.test_rho <- function(side="one",alpha=0.05,beta=0.2,delta=0.1){
  P <- switch(side,
              two=1-alpha/2,
              one=1-alpha)
  n <- ceiling((qnorm(P)+qnorm(1-beta))^2/delta^2)+3
  n
}

size_n.regII.test_rho_2 <- function(side="one",alpha=0.05,beta=0.2,rho0, delta){
  P <- switch(side,
              two=1-alpha/2,
              one=1-alpha)
  rho1 <- rho0-delta
  z0 <- 0.5*log((1+rho0)/(1-rho0))
  z1 <- 0.5*log((1+rho1)/(1-rho1))
  c <- z0-z1
  n <- ceiling((qnorm(P)+qnorm(1-beta))^2/c^2)+3
  n
}
