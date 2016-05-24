################################################################
##########          Densidades das SNI              ############

## Densidade/CDF da SN com locação escala #######

dSN <- function(y, mu , sigma2 = 1, shape=1){
  dens <- 2*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*((y - mu)/sqrt(sigma2)))
  return(dens)
}

dSNH <- function(rho, y, z, mu , sigma2 = 1, shape=1, rho.func){
  sigma2i<-nlwi(z,rho, rho.func)*sigma2
  dens <- 2*dnorm(y, mu, sqrt(sigma2i))*pnorm(shape*((y - mu)/sqrt(sigma2i)))
  return(dens)
}

## Densidade/CDF da ST com locação escala #######

dt.ls <- function(x, loc , sigma2 = 1,shape=1, nu = 4){
  d <- (x - loc)/sqrt(sigma2)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2)
  return(dens)
}

dtH.ls <- function(rho, y, z, loc , sigma2 = 1,shape=1, nu = 4, rho.func){

  sigma2i<-nlwi(z,rho, rho.func)*sigma2
  d <- (y - loc)/sqrt(sigma2i)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2i)
  return(dens)
}



## Densidade/CDF da Skew Normal Contaminada #######
  dSNC <- function(y, mu, sigma2, shape, nu){

    dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*sigma2^(-1/2)*(y-mu)))
    return(dens)
  }

    dSNCH <- function(rho, y, z, mu, sigma2, shape, nu, rho.func){
    sigma2i<-nlwi(z,rho, rho.func)*sigma2
    dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2i/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2i^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2i))*pnorm(shape*sigma2i^(-1/2)*(y-mu)))
    return(dens)
  }


### Densidade da Skew Slash  ######
dSS <- function(y, mu, sigma2, shape,nu){
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu[i],sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu[i]))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)

  }


dSSH <- function(rho, y, z, mu, sigma2, shape,nu, rho.func){
 sigma2i<-nlwi(z,rho, rho.func)*sigma2
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu[i],sqrt(sigma2i[i]/u))*pnorm(u^(1/2)*shape*(sigma2i[i]^(-1/2))*(y[i]-mu[i]))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}

##############################################################
######        FunÃ§ao  para o modelo Heterosced        #######

nlwi<-function(z,rho, type=1){
if(type == 1) resp<- exp(z*rho)
else resp<-z^rho
return(resp)
}
