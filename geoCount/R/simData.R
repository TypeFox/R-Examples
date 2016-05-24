####################################
#### Simulate locations
####################################
locGrid <- function(x, y, nx, ny){
  res <- as.matrix(expand.grid(seq(0, x, length=nx), seq(0, y, length=ny)))
  colnames(res) <- NULL
  res
}
locCircle <- function(r, np){
  agl <- seq(0, 2*pi, length=np+1)[1:np]
  cbind(r*cos(agl), r*sin(agl))
}
locSquad <- function(a, np){
  rbind(
  cbind(a, seq(-a, a, length=np)[1:(np-1)] ),
  cbind(seq(a, -a, length=np)[1:(np-1)], a ),
  cbind(-a, seq(a, -a, length=np)[1:(np-1)] ),
  cbind(seq(-a, a, length=np)[1:(np-1)], -a )
  )
}
####################################
#### Simulate data
####################################
simData <- function(loc, L=0, X=NULL, beta=0, cov.par, rho.family = "rhoPowerExp", Y.family="Poisson"){
  n <- nrow(loc)

  if(any(L==0)){
    L <- matrix(rep(1,n),,1)
    message("\nL contains zero!\nL is set to 1 for all locations!")
    } else { L <- matrix(L,,1)}
    
  U <- loc2U(loc)
  Z <- U2Z(U, cov.par, rho.family)
  D <- cbind( rep(1, n), X )
  mu.S <- D%*%beta
  z <- rnorm(n)
  S <- mu.S + chol(Z)%*%z
  if(Y.family=="Poisson"){
    Y <- rpois(n, L*exp(S))
  } else if(Y.family=="Binomial"){
          Y <- rbinom(n, L, exp(S)/(1+exp(S)))
          }
  list(data=Y, latent=S)
  }
####################################
#### END
####################################
