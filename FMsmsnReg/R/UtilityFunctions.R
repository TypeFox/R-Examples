###############################################################################
##########          Densidades das SNI    : paper 3          ################          #

## Densidade/CDF da SN com locacao escala #######
dSN <- function(y, mu, sigma2 = 1, shape=1)
{
  dens <- 2*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*((y - mu)/sqrt(sigma2)))
  return(dens)
}


## Densidade/CDF da ST com locacao escala #######
dt.ls <- function(x, loc, sigma2 = 1,shape=1, nu = 4)
{
  d    <- (x - loc)/sqrt(sigma2)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2)
  return(dens)
}


## Densidade/CDF da Skew Normal Contaminada #######
  dSNC <- function(y, mu, sigma2, shape, nu)
  {
   dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*sigma2^(-1/2)*(y-mu)))
   return(dens)
  }


### Densidade da Skew Slash  ######
 dSS <- function(y, mu, sigma2, shape,nu)
 {
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y))
  {
    f       <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu[i],sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu[i]))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}


###########    Densidades das Misturas de SNI   ##################
  d.mixedSN <- function(x, pi1, mu, sigma2, shape){
    # x: e o vetor de dados
    ## mu[,] uma matrix
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSN(x, mu[, j], sigma2[j], shape[j])
    return(dens)
  }

  d.mixedST <- function(x, pi1, mu, sigma2, shape, nu){
    # x: e o vetor de dados
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dt.ls(x, mu[,j], sigma2[j], shape[j], nu)
    return(dens)
  }

  d.mixedSNC <- function(x, pi1, mu, sigma2, shape, nu){
    # x: e o vetor de dados
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSNC(x, mu[,j], sigma2[j], shape[j], nu)
    return(dens)
  }

  d.mixedSS <- function(x, pi1, mu, sigma2, shape, nu){
    # x: e o vetor de dados
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSS(x, mu[,j], sigma2[j], shape[j], nu)
    return(dens)
  }

##########     FIM   Densidades das SNI             ############
################################################################

