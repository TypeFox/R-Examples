.simulate.rd <- function(n, known, unknown, N, mu=5, sigma=1, ...)
{
	#simulate d_i for n people using mu and sigma
	d <- rlnorm(n, mu, sigma)
	
	#create vector of p for binomial for simulate y_ik for each person i
	known.all <- c(known, unknown)
  p.pops <- known.all/N
	
	y <- matrix(nrow=n, ncol=length(known.all))
	for(i in 1:n)
	{
		di <- d[i]
		for(k in 1:length(known.all))
		{
			y[i,k] <- rbinom(1, size=as.integer(di), prob=p.pops[k])
		}
	}
	
	return(list(y=y, d=d))
}


#function for NK terms of log posterior (d_i integrated out)
# indices.k - takes index for one N_K
#NK - estimated value for N_K
#N - is population total
.loglik.NK.rd <- function(dat, indices.k, N, d, NK)
{	
	sum.di <- sum(d)
	sum.yik <- sum(dat[,indices.k])
	
	loglik <- sum.yik * log(NK/(N-NK)) + sum.di * log(1-NK/N) - log(NK)
	
	return(loglik)
}

#function for d_i log posterior terms
#takes into dat only for individual i

.loglik.di.rd <- function(dat, mu, sigma, di, N, known, NK)
{	
	sum.binom <- sum(lchoose(di, dat))
	
	sum.log.nk <- sum(log(1-known/N)) + sum(log(1-NK/N))
	
	loglik <- -log(di) - (log(di) - mu)^2 / (2*sigma^2) + sum.binom + di * sum.log.nk
	
	return(loglik)
}


#NKstart, NK.tuning, indices.k should all be same length (length of # of N_K being sampled)
#mu.prior, sigma.prior - length 2 vector having lower and upper limits of uniform prior

.mcmc.rd <- function(dat, known, N, indices.k, iterations, burnin, size, NK.start=killworth.start(dat,known,N)$NK.start, d.start=killworth.start(dat,known,N)$d.start, mu.start=killworth.start(dat,known,N)$mu.start, sigma.start=killworth.start(dat,known,N)$sigma.start, mu.prior=c(3,8), sigma.prior=c(1/4,2), NK.tuning=0.25*NK.start, d.tuning=0.25*d.start, ...)
{
  mu.values <- vector(length=size)
  sigma.values <- vector(length=size)
  NK.values <- matrix(nrow=length(indices.k), ncol=size)
  d.values <- matrix(nrow=dim(dat)[1], ncol=size)
  n <- dim(dat)[1]
  
  NK.curr <- NK.start
  d.curr <- d.start
  mu.curr <- mu.start
  sigma.curr <- sigma.start
  
  keep <- round(seq(from=burnin+1, to = burnin+iterations, length=size))
  keep.index <- 1
  
  for(t in 1:(burnin+iterations))
  {
    if(t %% 1000 == 0)
    {
      cat("Now on iteration ", t, "\n")
    }
    
    d <- d.curr
    NK <- NK.curr
    mu <- mu.curr
    sigma <- sigma.curr
    
    #############
    # compute d
    #############
    
    for(i in 1:n)
    {
      di.old <- d[i]
      di.new <- rnorm(1, mean=di.old, sd=d.tuning[i])
      
      if(di.new < max(dat[i,]))
      {
        d.curr[i] <- di.old
      } else
      {
        ratio <- .loglik.di.rd(dat=dat[i,], mu=mu, sigma=sigma, di=di.new, N=N, known=known, NK=NK) - .loglik.di.rd(dat=dat[i,], mu=mu, sigma=sigma, di=di.old, N=N, known=known, NK=NK)
        
        if(is.na(ratio))
        {
          cat("I broke for di for i = ", i, " where di.new was ", di.new, "\n")
        }
        
        accept <- min(0, ratio)
        
        logu <- log(runif(1, 0, 1))
        
        if(logu < accept)
        {
          d.curr[i] <- di.new
        } else
        {
          d.curr[i] <- di.old
        }
      }
    }
    
    d <- d.curr
    
    ###########
    ## update NK
    ###########
    
    for(j in 1:length(indices.k))
    {
      index.this.k <- indices.k[j]
      
      NK.old <- NK[j]
      NK.new <- rnorm(1, mean=NK.old, sd=NK.tuning[j])
      
      if(NK.new >= N || NK.new <= max(dat[, index.this.k]))
      {
        NK.curr[j] <- NK.old
      } else
      {
        ratio <- .loglik.NK.rd(dat, index.this.k, N, d, NK.new) - .loglik.NK.rd(dat, index.this.k, N, d, NK.old)
        
        #if(is.na(ratio))
        #{
        #	cat("I broke on NK where my K was ", indices.k, " and my new.NK was ", NK.new, "\n")
        #}
        
        accept <- min(0, ratio)
        
        logu <- log(runif(1, 0, 1))
        
        if(logu < accept)
        {
          NK.curr[j] <- NK.new
        } else
        {
          NK.curr[j] <- NK.old
        }
      }
    }
    
    NK <- NK.curr
    
    ###########
    ## update mu
    ###########
    
    mu.new <- rnorm(1, mean=sum(log(d))/n, sd=sigma/sqrt(n))
    
    # VERIFY  !!!!!!
    while(mu.new < mu.prior[1] || mu.new > mu.prior[2])
    {
      mu.new <- rnorm(1, mean=sum(log(d))/n, sd=sigma/sqrt(n))
    } 
    
    mu.curr <- mu.new
    
    mu <- mu.curr
    
    ###########
    ## update sigma
    ###########
    
    sigma.new <- sqrt(rinvgamma(1, shape=((n-1)/2), scale = 1/2*sum((log(d)-mu)^2)))
    
    while(sigma.new < sigma.prior[1] || sigma.new > sigma.prior[2])
    {
      sigma.new <- sqrt(rinvgamma(1, shape=((n-1)/2), scale = 1/2*sum((log(d)-mu)^2)))
    }
    
    sigma.curr <- sigma.new
    
    if(t == keep[keep.index])
    {
      mu.values[keep.index] <- mu.curr
      sigma.values[keep.index] <- sigma.curr
      NK.values[,keep.index] <- NK.curr
      d.values[,keep.index] <- d.curr
      keep.index <- keep.index+1
    }
    
  }
  
  return(list(NK.values = NK.values, d.values = d.values, mu.values = mu.values, sigma.values = sigma.values, iterations = iterations, burnin = burnin))
}