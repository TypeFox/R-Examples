#matrix of tau priors give a, b terms for beta distribution for tau
#each row is for a different subpopulation
#assumes unknown subpopulations that have tau's are last in "known"
.simulate.trans <- function(n, known, unknown, N, mu=5, sigma=1, tauK=rep(1, length(unknown)), ...)
{
	#simulate d_i for n people using mu and sigma
	d <- rlnorm(n, mu, sigma)
	
	#draw tau_K for unknown subpopulations
	tau <-c(rep(1, length(known)), tauK)
	
	#create vector of p for binomial for simulate y_ik for each person i
	known.all <- c(known, unknown)
	p.pops <- tau * known.all/N
	
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

#function for d_i log posterior terms
#takes into dat only for individual i
#tauK is tau_K for all K unknown
.loglik.di.trans <- function(dat, mu, sigma, di, N, known, wK)
{	
	sum.binom <- sum(lchoose(di, dat))
	
	sum.log.nk <- sum(log(1-known/N)) + sum(log(1-wK/N))
	
	loglik <- -log(di) - (log(di) - mu)^2 / (2*sigma^2) + sum.binom + di * sum.log.nk
	
	return(loglik)
}

#function for wK terms of log posterior where wK = NK * tauK
# indices.k - takes index for one N_K
# wK - estimated value for wK
# zK - estimated value for zK
# N - population total
# tauK.prior - prior terms for relevant tauK
.loglik.wK.trans <- function(dat, indices.k, N, d, wK, zK, tauK.prior)
{	
	aK <- tauK.prior[1]
	bK <- tauK.prior[2]
	
	sum.di <- sum(d)
	sum.yik <- sum(dat[,indices.k])
	
	loglik <- sum.yik * log(wK/(N-wK)) + sum.di * log(1-wK/N) + (aK - 2)/2 * log(wK) + (bK - 1)*log(1 - sqrt(wK/zK))
	
	return(loglik)
}

#function for zK terms of log posterior where zK = NK / tauK
# zK - estimated value for zK
# wK - estimated value for wK
# tauK.prior - prior terms for relevant tauK
.loglik.zK.trans <- function(wK, zK, tauK.prior)
{	
	aK <- tauK.prior[1]
	bK <- tauK.prior[2]
	
	loglik <- -(aK + 2)/2 * log(zK) + (bK - 1)*log(1 - sqrt(wK/zK))
	
	return(loglik)
}



#wKstart, wK.tuning, zKstart, zK.tuning, indices.k should all be same length (length of # of N_K being sampled)
#mu.prior, sigma.prior - length 2 vector having lower and upper limits of uniform prior

.mcmc.trans <- function(dat, known, N, indices.k, iterations, burnin, size, NK.start=killworth.start(dat,known,N)$NK.start, d.start=killworth.start(dat,known,N)$d.start, mu.start=killworth.start(dat,known,N)$mu.start, sigma.start=killworth.start(dat,known,N)$sigma.start, tauK.start=rep(0.5,length(indices.k)), mu.prior=c(3,8), sigma.prior=c(1/4,2), tauK.priors=matrix(rep(c(1,1),each=length(indices.k)),length(indices.k),2), wK.tuning=0.25*NK.start*tauK.start, zK.tuning=0.25*NK.start/tauK.start, d.tuning=0.25*d.start, ...)
{
	wK.start <- NK.start*tauK.start
  zK.start <- NK.start/tauK.start
  
  mu.values <- vector(length=size)
	sigma.values <- vector(length=size)
	d.values <- matrix(nrow=dim(dat)[1], ncol=size)
	wK.values <- matrix(nrow=length(indices.k), ncol=size)
	zK.values <- matrix(nrow=length(indices.k), ncol=size)
	n <- dim(dat)[1]
	
	d.curr <- d.start
	wK.curr <- wK.start
	zK.curr <- zK.start
	mu.curr <- mu.start
	sigma.curr <- sigma.start

  keep <- round(seq(from=burnin+1, to=burnin+iterations, length=size))
	keep.index <- 1
	
	for(t in 1:(burnin+iterations))
	{
		if(t %% 1000 == 0)
		{
			cat("Now on iteration ", t, "\n")
		}
		
		d <- d.curr
		wK <- wK.curr
		zK <- zK.curr
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
				ratio <- .loglik.di.trans(dat=dat[i,], mu=mu, sigma=sigma, di=di.new, N=N, known=known, wK=wK) - .loglik.di.trans(dat=dat[i,], mu=mu, sigma=sigma, di=di.old, N=N, known=known, wK=wK)
				
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
		## update wK
		###########
		
		for(j in 1:length(indices.k))
		{
			index.this.k <- indices.k[j]
			
			wK.old <- wK[j]
			wK.new <- rnorm(1, mean=wK.old, sd=wK.tuning[j])
			
			if(wK.new >= N || wK.new <= 0 || wK.new/zK[j] >= 1 || sqrt(wK.new * zK[j]) >= N) 
			{
				wK.curr[j] <- wK.old
			} else
			{	
				ratio <- .loglik.wK.trans(dat, index.this.k, N, d, wK.new, zK[j], tauK.priors[j,]) - .loglik.wK.trans(dat, index.this.k, N, d, wK.old, zK[j], tauK.priors[j,])
				
				if(is.na(ratio))
				{
					cat("I broke on wK where my K was ", indices.k, " and my new.wK was ", wK.new, " and zK is ", zK[j], "\n")
				}
			
				accept <- min(0, ratio)
			
				logu <- log(runif(1, 0, 1))
	  
	  			if(logu < accept)
	 			{
	   				wK.curr[j] <- wK.new
	  			} else
	 			{
	    			wK.curr[j] <- wK.old
	  			}
			}
		}
		
		wK <- wK.curr
		
		###########
		## update zK
		###########
		
		for(j in 1:length(indices.k))
		{			
			zK.old <- zK[j]
			zK.new <- rnorm(1, mean=zK.old, sd=zK.tuning[j])
						
			if(zK.new <= 0 || wK[j]/zK.new >= 1 || sqrt(wK[j] * zK.new) >= N) 
			{
				zK.curr[j] <- zK.old
			} else
			{			
				ratio <- .loglik.zK.trans(wK[j], zK.new, tauK.priors[j,]) - .loglik.zK.trans(wK[j], zK.old, tauK.priors[j,])
				
				if(is.na(ratio))
				{
					cat("I broke on zK where my K was ", index.this.k, " and my new.zK was ", zK.new, " and wK was currently ", wK[j], "\n")
				}
			
				accept <- min(0, ratio)
			
				logu <- log(runif(1, 0, 1))
	  
	  			if(logu < accept)
	 			{
	   				zK.curr[j] <- zK.new
	  			} else
	 			{
	    			zK.curr[j] <- zK.old
	  			}
	  		}
		}
	
		zK <- zK.curr
		
		###########
		## update mu
		###########

		mu.new <- rnorm(1, mean=sum(log(d))/n, sd=sigma/sqrt(n))
		
		while(mu.new < mu.prior[1] || mu.new > mu.prior[2])
		{
			cat("I'm in the mu while loops for mu.new of ", mu.new, "\n")
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
			#cat("I'm in the sigma while loops for sigma.new of ", sigma.new, "\n")
			sigma.new <- sqrt(rinvgamma(1, shape=((n-1)/2), scale = 1/2*sum((log(d)-mu)^2)))
		}
		
		sigma.curr <- sigma.new		
    
		if(t == keep[keep.index])
		{
		  d.values[,keep.index] <- d.curr
		  wK.values[,keep.index] <- wK.curr
		  zK.values[,keep.index] <- zK.curr
		  mu.values[keep.index] <- mu.curr
		  sigma.values[keep.index] <- sigma.curr
		  keep.index <- keep.index+1   
		}
	}
	
  NK.values <- sqrt(wK.values * zK.values)
  tauK.values <- sqrt(wK.values / zK.values)
	return(list(NK.values = NK.values, d.values = d.values, mu.values=mu.values, sigma.values=sigma.values, tauK.values = tauK.values, iterations = iterations, burnin = burnin))
}