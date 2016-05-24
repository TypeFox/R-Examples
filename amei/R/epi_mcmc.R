theta <- function(b, k, i)
{
  return(1 - (k/(k+b*i))^k)
}

mcmc.bknumu <- function(N, bkrate=1, b1, k1, nu1, mu1, itilde, rtilde, dtilde,
                      S, I, hyper)
{
  b <- k <- rep(0, N*bkrate+1)
  nu <- mu <- rep(0, N+1)
  b[1] <- b1; k[1] <- k1; nu[1] <- nu1; mu[1] <- mu1
  
  for(i in 1:(N*bkrate)){
    
    
    ## for each parameter, propose a new value, calculate
    ## the log posterior value of the current and proposed
    ## value (keeping the other parameters constant), and
    ## accept the proposed move with prob A (see below)  
    
    ## b
    u <- runif(1)
    q <- propose(b[i])
    b[i+1] <- q$b
    prob.old <- epi.log.post.b(b[i], itilde, k[i], S, I, hyper$bh)
    prob.new <- epi.log.post.b(b[i+1], itilde, k[i], S, I, hyper$bh)
    A <- exp( prob.new + q$lbak - prob.old - q$lfwd )
    if( u > A ) b[i+1] <- b[i] 
    
    ## k
    u <- runif(1)
    q <- propose(k[i])
    k[i+1] <- q$b
    prob.old <- epi.log.post.k(b[i+1], itilde, k[i], S, I, hyper$kh)
    prob.new <- epi.log.post.k(b[i+1], itilde, k[i+1], S, I, hyper$kh)
    A <- exp( prob.new + q$lbak - prob.old - q$lfwd )
    if( u > A ) k[i+1] <- k[i]
    
  }
  
  ## no longer doing gibbs steps or mcmc for mu, nu, we can sample directly
  nu <- epi.gibbs.nu(N+1, I, rtilde, dtilde, hyper$nuh)
  mu <- epi.gibbs.mu(N+1, I, rtilde, dtilde, hyper$muh)

  ## thin b & k by bkrate first
  b <- b[c(1,seq(2,(bkrate*N+1),length=N))]
  k <- k[c(1,seq(2,(bkrate*N+1),length=N))]

  ## then discard the first 10% as burn-in
  b <- b[(0.1*N):(N+1)]
  k <- k[(0.1*N):(N+1)]
  nu <- nu[(0.1*N):(N+1)]
  mu <- mu[(0.1*N):(N+1)]

  ## return the samples
  return(data.frame(b=b, k=k, nu=nu, mu=mu))
}

## used for a random walk proposal
propose <- function(b, l=3, h=4)
{
  b.new <- runif(1, l/h*b, h/l*b)
  fwd <- dunif(b.new, l/h*b, h/l*b)
  bak <- dunif(b, l/h*b.new, h/l*b.new)
  return(list(b=b.new, lfwd=log(fwd), lbak=log(bak)))
}	


## functions for calculating the appropriate posterior probs
epi.log.post.b <- function(b, itilde, k, S, I, ab)
{
  ## pretend there is at least one infected
  I[I<=0] <- 1

  ## calculate the binomial likelihood via theta, the
  ## probability of success
  t <- theta(b,k,I)
  llik <- sum(dbinom(itilde, S, t, log=TRUE))
  
  lprior <- dgamma(b, ab[1], scale=ab[2], log=TRUE) 

  ## check for problems in the likelihood
  if(!is.finite(llik) || !is.finite(lprior))
    stop("infinite likelihood or prior encountered")
        
  return( llik + lprior )	
}

epi.log.post.k <- function(b, itilde, k, S, I, ab)
{
  ## pretend there is at least one infected
  I[I<=0] <- 1

  ## calculate the binomial likelihood via theta, the
  ## probability of success
  t <- theta(b,k,I)
  llik <- sum(dbinom(itilde, S, t, log=TRUE))

  lprior <- dgamma(k, ab[1], scale=ab[2], log=TRUE) 

  ## check for problems in the likelihood
  if(!is.finite(llik) || !is.finite(lprior))
    stop("infinite likelihood or prior encountered")
  
  return( llik + lprior )	
}	


epi.gibbs.nu <- function(N, I, rtilde, dtilde, ab){

  alpha <- ab[1]
  beta <- ab[2]
  
  a2 <- alpha+sum(rtilde)
  b2 <- beta+sum(I)-sum(rtilde)
  pr <- rbeta(N, a2, shape2=b2)
  nu <- -log(1-pr)
  return(nu)
}

epi.gibbs.mu <- function(N, I, rtilde, dtilde, ab){

  alpha <- ab[1]
  beta <- ab[2]
  
  a2 <- alpha+sum(dtilde)
  b2 <- beta+sum(I)-sum(rtilde)-sum(dtilde)
  pd <- rbeta(N, a2, shape2=b2)
  mu <- -log(1-pd)
  return(mu)
}
