################################################################
##########          Densidades das SNI              ############

## Densidade/CDF da SN com locacao escala #######
dSN <- function(y, mu = 0, sigma2 = 1, shape=1){
  dens <- 2*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*((y - mu)/sqrt(sigma2)))
  return(dens)
}



## Densidade/CDF da ST com locacao escala #######
dt.ls <- function(x, loc = 0, sigma2 = 1,shape=1, nu = 4){
  d <- (x - loc)/sqrt(sigma2)
  dens <- 2*dt(d, df = nu)*pt(sqrt((1+nu)/(d^2+nu))*d*shape,1+nu)/sqrt(sigma2)
  return(dens)
}



## Densidade/CDF da Skew Normal Contaminada #######
  dSNC <- function(y, mu, sigma2, shape, nu){
    dens <- 2*(nu[1]*dnorm(y, mu, sqrt(sigma2/nu[2]))*pnorm(sqrt(nu[2])*shape*sigma2^(-1/2)*(y-mu)) + (1 - nu[1])*dnorm(y, mu, sqrt(sigma2))*pnorm(shape*sigma2^(-1/2)*(y-mu)))
    return(dens)
  }


### Densidade da Skew Slash  ######
dSS <- function(y, mu, sigma2, shape,nu){
  resp <- vector(mode = "numeric", length = length(y))
  for (i in 1:length(y)) {
    f <- function(u) 2*nu*u^(nu - 1)*dnorm(y[i],mu,sqrt(sigma2/u))*pnorm(u^(1/2)*shape*(sigma2^(-1/2))*(y[i]-mu))
    resp[i] <- integrate(f,0,1)$value
  }
  return(resp)
}


###########    Densidades das Misturas de SNI   ##################
  d.mixedSN <- function(x, pi1, mu, sigma2, shape){
    # x: eh o vetor de dados
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSN(x, mu[j], sigma2[j], shape[j])
    return(dens)
  }
  
  d.mixedST <- function(x, pi1, mu, sigma2, shape, nu){
    # x: eh o vetor de dados
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dt.ls(x, mu[j], sigma2[j], shape[j], nu)
    return(dens)
  }

  d.mixedSNC <- function(x, pi1, mu, sigma2, shape, nu){
    # x: eh o vetor de dados
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSNC(x, mu[j], sigma2[j], shape[j], nu)
    return(dens)
  }

  d.mixedSS <- function(x, pi1, mu, sigma2, shape, nu){
    # x: eh o vetor de dados
    # outros parametros devem ser do tipo vetor c() de dimensao g (qtd de misturas)
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dSS(x, mu[j], sigma2[j], shape[j], nu)
    return(dens)
  }

##########     FIM   Densidades das SNI             ############
################################################################


#######################################################################
########      Algoritmo EM e funcoes MULTIVARIADAS         ############
#######################################################################

#______________________________________________________________________________________________
################################################################
##########          Densidades das SNI              ############
#pacote necessario
##require(mvtnorm)
##require(mnormt)

#matrix.sqrt <- function(A) {
#  sva <- svd(A)
#  if (min(sva$d)>=0)
#    Asqrt <- t(sva$v %*% (t(sva$u) * sqrt(sva$d)))
#  else
#    stop("Matrix square root is not defined")
#  return(Asqrt)
#}


## Densidade/CDF da SN com locacao escala #######
dmvSN <- function(y, mu, Sigma, lambda){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensao ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimensao igual a ncol(y) = p
  #Sigma: Matrix p x p
  n <- nrow(y)
  p <- ncol(y)
  dens <- 2*dmvnorm(y, mu, Sigma)*pnorm( apply( matrix(rep(t(lambda)%*%solve(matrix.sqrt(Sigma)),n), n, p, byrow = TRUE)*(y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1, sum  ) )
  return(dens)
}


## Densidade/CDF da ST com locacao escala #######
dmvt.ls <- function(y, mu, Sigma, lambda, nu){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensao ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimensao igual a ncol(y) = p
  #Sigma: Matrix p x p
  n <- nrow(y)
  p <- ncol(y)
  denst <- (gamma((p+nu)/2)/(gamma(nu/2)*pi^(p/2)))*nu^(-p/2)*det(Sigma)^(-1/2)*(1 + mahalanobis(y, mu, Sigma)/nu)^(-(p+nu)/2)
  dens <- 2*(denst)*pt(sqrt((p + nu)/(mahalanobis(y, mu, Sigma) + nu))*apply( matrix(rep(t(lambda)%*%solve(matrix.sqrt(Sigma)),n), n, p, byrow = TRUE)*(y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1, sum  )    , df = nu + p   )
  return(dens)
}


## Densidade/CDF da Skew Normal Contaminada #######
dmvSNC <- function(y, mu, Sigma, lambda, nu){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensao ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimensao igual a ncol(y) = p
  #Sigma: Matrix p x p
  n <- nrow(y)
  p <- ncol(y)
  dens <- 2*(nu[1]*dmvnorm(y, mu, Sigma/nu[2])*pnorm(sqrt(nu[2])*apply( matrix(rep(t(lambda)%*%solve(matrix.sqrt(Sigma)),n), n, p, byrow = TRUE)*(y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1, sum  ) ) + (1 - nu[1])*dmvnorm(y, mu, Sigma)*pnorm(apply( matrix(rep(t(lambda)%*%solve(matrix.sqrt(Sigma)),n), n, p, byrow = TRUE)*(y - matrix(rep(mu, n), n, p, byrow = TRUE)), 1, sum  )) )
  return(dens)
}


dmvSS <- function(y, mu, Sigma, lambda, nu){
  #y: deve ser uma matrix onde cada linha tem um vetor de dados multivariados de dimensao ncol(y) = p. nrow(y) = tamanho da amostra
  #mu, lambda: devem ser do tipo vetor de mesma dimensao igual a ncol(y) = p
  #Sigma: Matrix p x p
  if (is.matrix(y)) {
    n <- nrow(y)
    p <- ncol(y)
    resp <- c()
    for (i in 1:n) {
      di <- mahalanobis(y[i,], mu, Sigma)
      f <- function(u) 2*nu*u^(nu - 1)*( (u/(2*pi))^(p/2)*det(Sigma)^(-1/2)*exp(-u*di/2))*pnorm(u^(1/2)*t(lambda)%*%solve(matrix.sqrt(Sigma))%*%as.matrix(y[i,] - mu) )
      resp[i] <- integrate(f,0,1)$value
    }

  } else {
    n <- 1
    p <- length(mu)
    resp <- c()
    y <- matrix(y, 1, p)
    for (i in 1:n) {
      di <- mahalanobis(y[i,], mu, Sigma)
      f <- function(u) 2*nu*u^(nu - 1)*( (u/(2*pi))^(p/2)*det(Sigma)^(-1/2)*exp(-u*di/2))*pnorm(u^(1/2)*t(lambda)%*%solve(matrix.sqrt(Sigma))%*%(y[i,] - mu) )
      resp[i] <- integrate(f,0,1)$value
    }
  }
  return(resp)
}


###########    Densidades das Misturas de SNI   ##################
  d.mixedmvSN <- function(y, pi1, mu, Sigma, lambda){
    #y: eh a matriz de dados
    #pi1: deve ser do tipo vetor de dimensao g
    #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
    #lambda: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dmvSN(y, mu[[j]], Sigma[[j]], lambda[[j]])
    return(dens)
  }

  d.mixedmvST <- function(y, pi1, mu, Sigma, lambda, nu){
    # y: eh a matriz de dados
    #pi1: deve ser do tipo vetor de dimensao g
    #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
    #lambda: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    #nu: eh um numero
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dmvt.ls(y, mu[[j]], Sigma[[j]], lambda[[j]], nu)
    return(dens)
  }

  d.mixedmvSNC <- function(y, pi1, mu, Sigma, lambda, nu){
    # y: eh a matriz de dados
    #pi1: deve ser do tipo vetor de dimensao g
    #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
    #lambda: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    #nu: eh um vetor de tamanho 2
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dmvSNC(y, mu[[j]], Sigma[[j]], lambda[[j]], nu)
    return(dens)
  }


  d.mixedmvSS <- function(y, pi1, mu, Sigma, lambda, nu){
    # y: eh a matriz de dados
    #pi1: deve ser do tipo vetor de dimensao g
    #mu: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    #Sigma: deve ser do tipo list com g entradas. Cada entrada do list deve ser uma matriz p x p
    #lambda: deve ser do tipo list com g entradas. Cada entrada do list deve ser um vetor de dimensao p
    #nu: eh um vetor de tamanho 2
    g <- length(pi1)
    dens <- 0
    for (j in 1:g) dens <- dens + pi1[j]*dmvSS(y, mu[[j]], Sigma[[j]], lambda[[j]], nu)
    return(dens)
  }


##########        FIM  - Densidades das SNI         ############
################################################################
