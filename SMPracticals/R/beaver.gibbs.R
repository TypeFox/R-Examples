"beaver.gibbs" <-
function(init, y, R=10, a=1, b=0.05)
{  # Gibbs sampler starting from init with R iterations
  sim.g <- function(y, mu1, mu2, lam)
  {  
    n <- length(y)
    p <- rep(NA,n-1)
    for (g in 1:(n-1)) 
     { mu <- c(rep(mu1,g),rep(mu2,n-g))
       p[g] <- sum((y-mu)^2) }
    sample(1:(n-1), size=1, prob=exp(-lam*p/2))  
  }
  lik <- function(y, para, a, b)
  {  # computes log likelihood, log prior at current parameters
    mu1 <- para[1]
    mu2 <- para[2]
    lam <- para[3]
    g <- para[4]
    n <- length(y)
    mu <- c(rep(mu1,g),rep(mu2,n-g))
    L <- sum( dnorm(y, mu, sd=sqrt(1/lam), log=TRUE) )  
    P <- dgamma(lam, shape=a, rate=b, log=TRUE)         
    c(L, P)
  }
  step <- function(j, para, y, a, b)
  {  # single step of random scan Gibbs sampler: update para[j]
    mu1 <- para[1]
    mu2 <- para[2]
    lam <- para[3]
    g <- para[4]
    n <- length(y)
    mu <- c(rep(mu1,g),rep(mu2,n-g))
    SS <- sum((y-mu)^2)
    sim <- switch(j, 
             rnorm(1, mean(y[1:g]), 1/sqrt(lam*g)),           # mu1
             rnorm(1, mean(y[(g+1):n]), 1/sqrt(lam*(n-g))),   # mu2
             rgamma(1, shape=a+n/2, rate=b+SS/2),             # lambda
             sim.g(y, mu1, mu2, lam))                         # gamma
    para[j] <- sim
    para
  }

 paras.out <- matrix(init, R, length(init), byrow=TRUE)
  log.lik <- lik(y, init, a, b)
  lik.out <- matrix(log.lik, R, length(log.lik), byrow=TRUE)
  for (r in 2:R)
  {  j <- sample(1:4)  # random order for four possible updates
     para <- paras.out[r-1,]
     for (i in 1:length(j)) para <- step(j[i], para, y, a, b)
     paras.out[r,] <- para
     lik.out[r,] <- lik(y, para, a, b)  }
  cbind(paras.out, lik.out)   # output is parameters, log likelihood, log prior
}

