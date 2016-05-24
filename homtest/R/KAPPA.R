# 2005-09-23, Alberto Viglione
#
# KAPPA provides the link between L-moments of a sample and the four parameter
# kappa distribution.
# Algorithm based on Hosking and Wallis...

# h=-1 => is the generalized logistic distribution
# h=0 => is the generalized extreme value distribution
# h=1 => is the generalized Pareto distribution


f.kappa <- function (x,xi,alfa,k,h) {

  if (k > 0) {
    if (x > (xi + alfa/k)) {
      stop("if k>0 x must be lower than xi + alfa/k")
    } 
  }
  
  if (h > 0) {
    if (x < (xi + alfa*(1 - h^(-k))/k)) {
      stop("if h>0 x must be higher than xi + alfa*(1 - h^(-k))/k")
    } 
  }
  else if (k < 0) {
    if (x < (xi + alfa/k)) {
      stop("if h<=0 and k<0 x must be higher than xi + alfa/k")
    } 
  }

  if (k == 0) {
    k <- 10^(-100)
  } 

  if (h == 0) {
    f <- f.GEV(x,xi,alfa,k)
  } 
  else {
    f <- alfa^(-1) *(1 - k*(x - xi)/alfa)^(1/(k-1)) * (F.kappa(x,xi,alfa,k,h))^(1/h)
  }

  return(f)
}

F.kappa <- function (x,xi,alfa,k,h) {

  if (k > 0) {
    if (x > (xi + alfa/k)) {
      stop("if k>0 x must be lower than xi + alfa/k")
    } 
  }
  
  if (h > 0) {
    if (x < (xi + alfa*(1 - h^(-k))/k)) {
      stop("if h>0 x must be higher than xi + alfa*(1 - h^(-k))/k")
    } 
  }
  else if (k < 0) {
    if (x < (xi + alfa/k)) {
      stop("if h<=0 and k<0 x must be higher than xi + alfa/k")
    } 
  }

  if (k == 0) {
    k <- 10^(-100)
  } 

  if (h == 0) {
    F <- F.GEV(x,xi,alfa,k)
  } 
  else {
    F <- (1 - h*(1 - k*(x - xi)/alfa)^(1/k))^(1/h)
  }
  
  return(F)
}

invF.kappa <- function (F,xi,alfa,k,h) {

  if ((F < 0) || (F > 1)) {
    stop("F must be between 0 and 1")
  } 


  if (k == 0) {
    k <- 10^(-100)
  } 
    
  if (h == 0) {
    x <- invF.GEV(F,xi,alfa,k)
  } 
  else {
    x <- xi + (alfa/k) * (1 - ((1 - F^h)/h)^k)
  }

  return(x)
}

Lmom.kappa <- function(xi,alfa,k,h) {

  xi <- as.numeric(xi)
  alfa <- as.numeric(alfa)
  k <- as.numeric(k)
  h <- as.numeric(h)
  
  quanti <- length(k)
  lambda1 <- rep(NA,quanti)
  lambda2 <- rep(NA,quanti)
  tau3 <- rep(NA,quanti)
  tau4 <- rep(NA,quanti)
  for (i in 1:quanti) {
    if (((k[i] < -1) && (h[i] >= 0)) || ((h[i] < 0) && ((k[i] <= -1) || (k[i] >= -1/h[i])))) {
      stop("L-moments are defined if h>=0 and k>-1, or if h<0 and -1<k<-1/h")
    } 

    if (k[i] == 0) {
      k[i] <- 10^(-100)
    } 

    if (h[i] == 0) {
      output <- Lmom.GEV(xi[i],alfa[i],k[i])
    } 
    else {
      g <- rep(NA,4)
      for (r in 1:4) {
        if (h[i] > 0) {
          g[r] <- (r*gamma(1+k[i])*gamma(r/h[i])) / (h[i]^(1+k[i]) *gamma(1+k[i]+r/h[i]))
        }
        else {
          g[r]=(r*gamma(1+k[i])*gamma(-k[i]-r/h[i])) / ((-h[i])^(1+k[i]) *gamma(1-r/h[i]))
        }
      }
  
      lambda1[i] <- xi[i] + alfa[i]*(1 - g[1])/k[i] 
      lambda2[i] <- alfa[i]*(g[1] - g[2])/k[i] 
      tau3[i] <- (-g[1] + 3*g[2] -2*g[3])/(g[1]-g[2])
      tau4[i] <- -(-g[1] + 6*g[2] -10*g[3] + 5*g[4])/(g[1]-g[2])
    }
  }
  output <- list(lambda1=lambda1, lambda2=lambda2, tau3=tau3, tau4=tau4)

  return(output)

}

par.kappa <- function(lambda1,lambda2,tau3,tau4) {
  
  lambda1 <- as.numeric(lambda1)
  lambda2 <- as.numeric(lambda2)
  tau3 <- as.numeric(tau3)
  tau4 <- as.numeric(tau4)

  sumquad.tau3tau4 = function (k.h,t3.t4) {

    k <- k.h[1]
    h <- k.h[2]
    t3 <- t3.t4[1]
    t4 <- t3.t4[2]
    
    error <- FALSE
    if (((k < -1) && (h >= 0)) || ((h < 0) && ((k <= -1) || (k >= -1/h)))) {
      stop("L-moments are defined if h>=0 and k>-1, or if h<0 and -1<k<-1/h")
      error=TRUE
    }

    g <- c(0,0,0,0)
  
    if (h == 0) {
      # GEV
      tau3 <- 2*(1 - 3^(-k))/(1 - 2^(-k)) - 3
      tau4 <- (5*(1 - 4^(-k)) - 10*(1 - 3^(-k)) + 6*(1 - 2^(-k)))/(1 - 2^(-k))
    }
    else {
      for (r in 1:4) {
        if (h > 0) {
          g[r] <- (r*gamma(1+k)*gamma(r/h)) / (h^(1+k) *gamma(1+k+r/h))
        }
        else {
          g[r]=(r*gamma(1+k)*gamma(-k-r/h)) / ((-h)^(1+k) *gamma(1-r/h))
        }
      }
  
      tau3 <- (-g[1] + 3*g[2] -2*g[3])/(g[1]-g[2])
      tau4 <- -(-g[1] + 6*g[2] -10*g[3] + 5*g[4])/(g[1]-g[2])
    }
    
    if (error == FALSE) { 
      output <- (t3-tau3)^2 + (t4-tau4)^2
    }
    else { 
      output <- -1
    }
  
    return(output)  
  }

  xi.alfa = function (lambda1,lambda2,k,h) {

    if (((k < -1) && (h >= 0)) || ((h < 0) && ((k <= -1) || (k >= -1/h)))) {
      stop("L-moments are defined if h>=0 and k>-1, or if h<0 and -1<k<-1/h")
    }

    g <- c(0,0)
    if (h == 0) {
      # GEV
      alfa <- (lambda2*k)/((1 - 2^(-k))*gamma(1+k))
      xi <- lambda1 - alfa*(1 - gamma(1+k))/k
    }
    else {
      for (r in 1:2) {
        if (h > 0) {
          g[r] <- (r*gamma(1+k)*gamma(r/h)) / (h^(1+k) *gamma(1+k+r/h))
        }
        else {
          g[r]=(r*gamma(1+k)*gamma(-k-r/h)) / ((-h)^(1+k) *gamma(1-r/h))
        }
      }
  
      alfa <- (lambda2*k)/(g[1]-g[2])
      xi <- lambda1 - alfa*(1-g[1])/k
    }

    output <- list(xi=xi, alfa=alfa)

    return(output)  
  }

  quanti <- length(tau3)
  k <- rep(NA,quanti)
  h <- rep(NA,quanti)
  xi <- rep(NA,quanti)
  alfa <- rep(NA,quanti)
  for (i in 1:quanti) {
    minimo <- optim(c(1,1),sumquad.tau3tau4,t3.t4=c(tau3[i],tau4[i]))
	if (minimo$value != -1) {
	  k[i] <- minimo$par[1]
	  h[i] <- minimo$par[2]
	  pp <- xi.alfa(lambda1[i],lambda2[i],k[i],h[i])
	  xi[i] <- pp$xi
	  alfa[i] <- pp$alfa
	}
  }
  output <- list(xi=xi,alfa=alfa,k=k,h=h)
  return(output)
}

rand.kappa <- function(numerosita,xi,alfa,k,h) {

  F <- runif(numerosita, min=0.0000000001, max=0.9999999999)
  x <- invF.kappa(F,xi,alfa,k,h)

  return(x)
}

