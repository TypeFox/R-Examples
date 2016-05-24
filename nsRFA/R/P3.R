# 2005-09-23, Alberto Viglione
#
# A.9) Pearson type III distribution

# gamm=0 => normal, gamm=2 => exponential, gamm=-2 => reverse exponential distribution

f.gamma <- function (x,xi,beta,alfa) {

  #f <- ((x - xi)^(alfa - 1) * exp(-(x - xi)/beta))/(abs(beta)*beta^(alfa-1) * gamma(alfa))

  if (beta > 0) {
    f <- dgamma(x-xi,shape=alfa,scale=beta)
  }
  else {
    f <- dgamma(xi-x,shape=alfa,scale=-beta)
  }

 return(f)
}

F.gamma <- function (x,xi,beta,alfa) {

  if (beta > 0) {
    #F <- pgamma((x - xi)/beta, alfa)
    F <- pgamma(x-xi,shape=alfa,scale=beta)
  }
  else {
    #F <- 1 - pgamma((xi - x)/beta, alfa) 
    F <- 1 - pgamma(xi-x,shape=alfa,scale=-beta)
  }

  return(F)
}

invF.gamma <- function (F,xi,beta,alfa) {

  # if ((F < 0) || (F > 1)) {
  #   stop("F must be between 0 and 1")
  # } 

  if (beta >= 0) {
    x.st <- qgamma(F, shape=alfa,scale=beta)
    x <- x.st + xi
  }
  else {
    x.st <- qgamma(1-F, shape=alfa,scale=-beta)
    x <- xi - x.st
  }

  return(x)
}

Lmom.gamma <- function(xi,beta,alfa) {

  A0 = 0.32573501
  A1 = 0.16869150
  A2 = 0.078327243
  A3 = -0.0029120539
  B1 = 0.46697102
  B2 = 0.24255406
  C0 = 0.12260172
  C1 = 0.053730130
  C2 = 0.043384378
  C3 = 0.011101277
  D1 = 0.18324466
  D2 = 0.20166036
  E1 = 2.3807576
  E2 = 1.5931792
  E3 = 0.11618371
  F1 = 5.1533299
  F2 = 7.1425260
  F3 = 1.9745056
  G1 = 2.1235833
  G2 = 4.1670213
  G3 = 3.1925299
  H1 = 9.0551443
  H2 = 26.649995
  H3 = 26.193668
  
  quanti <- length(beta)
  lambda1 <- rep(NA,quanti)
  lambda2 <- rep(NA,quanti)
  tau3 <- rep(NA,quanti)
  tau4 <- rep(NA,quanti)
  for (i in 1:quanti) {
    if (beta[i] >= 0) {
      lambda1[i] <- xi[i] + alfa[i]*beta[i]
      lambda2[i] <- abs(pi^(-0.5) *beta[i] * gamma(alfa[i] + 0.5)/gamma(alfa[i]))
      # tau3 <- 6 * pbeta(1/3,alfa,2*alfa) - 3
      if (alfa[i] >= 1) {
        tau3[i] <- alfa[i]^(-0.5) * (A0 + A1*alfa[i]^(-1) + A2*alfa[i]^(-2) + A3*alfa[i]^(-3))/
                   (1 + B1*alfa[i]^(-1) + B2*alfa[i]^(-2))
        tau4[i] <- (C0 + C1*alfa[i]^(-1) + C2*alfa[i]^(-2) + C3*alfa[i]^(-3))/(1 + D1*alfa[i]^(-1) + D2*alfa[i]^(-2))
      }
      else {
        tau3[i] <- (1 + E1*alfa[i] + E2*alfa[i]^2 + E3*alfa[i]^3)/(1 + F1*alfa[i] + F2*alfa[i]^2 + F3*alfa[i]^3)
        tau4[i] <- (1 + G1*alfa[i] + G2*alfa[i]^2 + G3*alfa[i]^3)/(1 + H1*alfa[i] + H2*alfa[i]^2 + H3*alfa[i]^3)
      }
    }
    else {
      lambda1[i] <- xi[i] + alfa[i]*beta[i]
      lambda2[i] <- abs(pi^(-0.5) *beta[i] * gamma(alfa[i] + 0.5)/gamma(alfa[i]))
      if (alfa[i] >= 1) {
        tau3[i] <- -(alfa[i]^(-0.5) * (A0 + A1*alfa[i]^(-1) + A2*alfa[i]^(-2) + A3*alfa[i]^(-3))/
                   (1 + B1*alfa[i]^(-1) + B2*alfa[i]^(-2)))
        tau4[i] <- (C0 + C1*alfa[i]^(-1) + C2*alfa[i]^(-2) + C3*alfa[i]^(-3))/(1 + D1*alfa[i]^(-1) + D2*alfa[i]^(-2))
      }
      else {
        tau3[i] <- -(1 + E1*alfa[i] + E2*alfa[i]^2 + E3*alfa[i]^3)/(1 + F1*alfa[i] + F2*alfa[i]^2 + F3*alfa[i]^3)
        tau4[i] <- (1 + G1*alfa[i] + G2*alfa[i]^2 + G3*alfa[i]^3)/(1 + H1*alfa[i] + H2*alfa[i]^2 + H3*alfa[i]^3)
      }
    }
  }
  
  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)
  
  return(output)
}

#par.gamma <- function(lambda1,lambda2,tau3) {
#
#  lambda1 <- as.numeric(lambda1)
#  lambda2 <- as.numeric(lambda2)
#  tau3 <- as.numeric(tau3)
#
#  if ((abs(tau3) > 0)&&(abs(tau3) < 1/3)) {
#   z <- 3*pi*tau3^2
#   alfa <- (1 + 0.2906*z)/(z + 0.1882*z^2 + 0.0442*z^3)
#  }
#  else if ((abs(tau3) >= 1/3)&&(abs(tau3) < 1)) {
#   z <- 1 - abs(tau3)
#   alfa <- (0.36067*z - 0.59567*z^2 + 0.25361*z^3)/(1 - 2.78861*z + 2.56096*z^2 - 0.77045*z^3)
#  }
#  if(alfa<100) {
#   sigma <- lambda2*pi^(0.5) * alfa^(0.5) * gamma(alfa)/gamma(alfa + 0.5)
#   beta <- 0.5*sigma*abs(2*alfa^(-0.5))
#   beta <- sign(tau3)*beta
#   xi <- lambda1 - alfa*beta
#   output <- list(xi=xi, beta=beta, alfa=alfa)
#  }
#  else {
#   mu <- lambda1
#   sigma <- sqrt(pi)*lambda2/(1-1/(8*alfa)+1/(128*alfa^2))
#   output <- list(mu=mu, sigma=sigma)
#  }
#
#  return(output)
#}

par.gamma <- function(lambda1,lambda2,tau3) {

  lambda1 <- as.numeric(lambda1)
  lambda2 <- as.numeric(lambda2)
  tau3 <- as.numeric(tau3)
  quanti <- length(tau3)

  xis <- rep(NA, quanti)
  betas <- rep(NA, quanti)
  alfas <- rep(NA, quanti)
  mus <- rep(NA, quanti)
  sigmas <- rep(NA, quanti)
  for (i in 1:quanti) {
   if ((abs(tau3[i]) > 0)&&(abs(tau3[i]) < 1/3)) {
    z <- 3*pi*tau3[i]^2
    alfa <- (1 + 0.2906*z)/(z + 0.1882*z^2 + 0.0442*z^3)
   }
   else if ((abs(tau3[i]) >= 1/3)&&(abs(tau3[i]) < 1)) {
    z <- 1 - abs(tau3[i])
    alfa <- (0.36067*z - 0.59567*z^2 + 0.25361*z^3)/(1 - 2.78861*z + 2.56096*z^2 - 0.77045*z^3)
   }
   if(alfa<100) {
    sigma <- lambda2[i]*pi^(0.5) * alfa^(0.5) * gamma(alfa)/gamma(alfa + 0.5)
    beta <- 0.5*sigma*abs(2*alfa^(-0.5))
    beta <- sign(tau3[i])*beta
    xi <- lambda1[i] - alfa*beta
    xis[i] <- xi; betas[i] <- beta; alfas[i] <- alfa
   }
   else {
    mu <- lambda1[i]
    sigma <- sqrt(pi)*lambda2[i]/(1-1/(8*alfa)+1/(128*alfa^2))
    mus[i] <- mu; sigmas[i] <- sigma
   }
  }
  output <- list(xi=xis, beta=betas, alfa=alfas, mu=mus, sigma=sigmas)
  return(output)
}


rand.gamma <- function(numerosita,xi,beta,alfa) {

  #F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  #x <- invF.gamma(F,xi,beta,alfa)
  #x <- xi + beta*rgamma(numerosita, shape=alfa)
  if (beta>0) {
   x <- xi + rgamma(numerosita, shape=alfa, scale=beta)
  }
  else {
   x <- xi - rgamma(numerosita, shape=alfa, scale=-beta)
  }
  return(x)
}

mom2par.gamma <- function(mu,sigma,gamm) {
  
  if(gamm==0) {stop("The distribution is Normal")}
  else {
   alfa <- 4/(gamm^2)
   beta <- 0.5*sigma*abs(gamm) 
   xi <- mu - 2*sigma/gamm
  }
  
  output <- list(alfa=alfa, beta=beta, xi=xi)

  return(output)
}

par2mom.gamma <- function(alfa,beta,xi) {

  gamm <- 2/sqrt(alfa)
  sigma <- 2*beta/abs(gamm)
  mu <- xi + 2*sigma/gamm

  output <- list(mu=mu, sigma=sigma, gamm=gamm)

  return(output)
}

