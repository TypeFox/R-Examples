"hus.gibbs" <-
function(init, y, R=10, a1=1, a2=1, c=0.01, d=0.01)
{  # Gibbs sampler starting from init with R iterations
  sim.tau <- function(y, lam1, lam2)
  {  
    s <- cumsum(y)
    tau <- 1:length(y)
    p <- exp( s*log(lam1/lam2) + tau*(lam2-lam1) )
    sample(1:length(y), size=1, prob=p)  
  }
  lik <- function(y, para, a1, a2, c, d)
  {  # computes log likelihood, log prior at current parameters
    lam1 <- para[1]
    lam2 <- para[2]
    beta1 <- para[3]
    beta2 <- para[4]
    tau <- para[5]
    lam <- c(rep(lam1,tau),rep(lam2,length(y)-tau))
    L <- sum( dpois(y, lam, log=TRUE) )  
    P <- dgamma(lam1, shape=a1, rate=beta1, log=TRUE) +  
         dgamma(lam2, shape=a2, rate=beta2, log=TRUE) +
         dgamma(beta1, shape=c, rate=d, log=TRUE) +
         dgamma(beta2, shape=c, rate=d, log=TRUE)
    c(L, P)
  }
  step <- function(j, para, y, a1, a2, c, d)
  {  # single step of random scan Gibbs sampler: update para[j]
    lam1 <- para[1]
    lam2 <- para[2]
    beta1 <- para[3]
    beta2 <- para[4]
    tau <- para[5]
    s.tau <- sum(y[1:tau])
    s.n <- sum(y)
    n <- length(y)
    sim <- switch(j, 
             rgamma(1, shape=a1+s.tau, rate=tau+beta1),            # lam1
             rgamma(1, shape=a2+s.n-s.tau, rate=n-tau+beta2),      # lam2
             rgamma(1, shape=a1+c, rate=lam1+d),                   # beta1
             rgamma(1, shape=a2+c, rate=lam2+d),                   # beta2
             sim.tau(y, lam1, lam2))                               # tau
    para[j] <- sim
    para
  }
  paras.out <- matrix(init, R, length(init), byrow=TRUE)
  log.lik <- lik(y, init, a1, a2, c, d)
  lik.out <- matrix(log.lik, R, length(log.lik), byrow=TRUE)
  for (r in 2:R)
  {  j <- sample(1:5)  # random order for five possible updates
     para <- paras.out[r-1,]
     for (i in 1:length(j)) para <- step(j[i], para, y, a1, a2, c, d)
     paras.out[r,] <- para
     lik.out[r,] <- lik(y, para, a1, a2, c, d)  }
  cbind(paras.out, lik.out)   # output is parameters, log likelihood, log prior
}

