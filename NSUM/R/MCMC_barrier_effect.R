.simulate.bar <- function(n, known, unknown, N, mu=5, sigma=1, rho=rep(0.1, length(c(known, unknown))), ...)
{
	#simulate d_i for n people using mu and sigma
	d <- rlnorm(n, mu, sigma)
	known.all <- c(known, unknown)
  
	#simulate per person proportions using beta distribution
	proportions <- matrix(nrow=n, ncol=length(known.all))
	for(k in 1:length(known.all))
	{
		a <- (known.all[k]/N)*(1/rho[k] - 1)
		b <- (1-known.all[k]/N)*(1/rho[k] - 1)
		
		proportions[,k]	<- rbeta(n, shape1=a, shape2=b)
	}
	
	y <- matrix(nrow=n, ncol=length(known.all))
	for(i in 1:n)
	{
		di <- d[i]
		for(k in 1:length(known.all))
		{
			y[i,k] <- rbinom(1, size=as.integer(di), prob=proportions[i,k])
		}
	}
	
	return(list(y=y, d=d))
}


#function for mK terms of log posterior 
# indices.k - takes index for one m_K
#mK - estimated value for m_K
#n - is number of respondents in survey
#rhoK - rho for a given mK
#index - k of term being updated
.loglik.rhok.bar <- function(dat, index, n, d, mk, rhok)
{
	sum.lbeta1 <- sum(lbeta(mk*(1/rhok - 1) + dat[,index], d + (1-mk)*(1/rhok - 1) - dat[,index]))
	
	sum.lbeta2 <- n*lbeta(mk*(1/rhok - 1), (1-mk)*(1/rhok - 1))
	
	loglik <- sum.lbeta1 - sum.lbeta2 
		
	return(loglik)
}




#function for rho_k terms of log posterior 
#rhok - estimated value for rho_k
#mk - single m_k value (rho_k matches for the k)
#index.k - index of K being updated
.loglik.mK.bar <- function(dat, index.k, n, d, mK, rhoK)
{
	sum.lbeta1 <- sum(lbeta(mK*(1/rhoK - 1) + dat[,index.k], d + (1-mK)*(1/rhoK - 1) - dat[,index.k]))
	
	sum.lbeta2 <- n*lbeta(mK*(1/rhoK - 1), (1-mK)*(1/rhoK - 1))
	
	loglik <- sum.lbeta1 - sum.lbeta2 - log(mK)
	
	return(loglik)
}


#function for d_i log posterior terms
#takes into dat only for individual i

.loglik.di.bar <- function(dat, mu, sigma, di, mknown, mK, rho)
{	
	m.all <- c(mknown, mK)
	
	sum.binom <- sum(lchoose(di, dat))
	
	sum.lgamma1 <- sum(lgamma(di+(1-m.all)*(1/rho - 1) - dat))
	sum.lgamma2 <- sum(lgamma(di+(1/rho - 1)))
	
	loglik <- -log(di) - (log(di) - mu)^2 / (2*sigma^2) + sum.binom + sum.lgamma1 - sum.lgamma2
	
	return(loglik)
}

#NKstart, NK.tuning, indices.k should all be same length (length of # of N_K being sampled)
#mu.prior, sigma.prior - length 2 vector having lower and upper limits of uniform prior

.mcmc.bar <- function(dat, known, N, indices.k, iterations, burnin, size, NK.start=killworth.start(dat,known,N)$NK.start, d.start=killworth.start(dat,known,N)$d.start, mu.start=killworth.start(dat,known,N)$mu.start, sigma.start=killworth.start(dat,known,N)$sigma.start, rho.start=rep(0.1,dim(dat)[2]), mu.prior=c(3,8), sigma.prior=c(1/4,2), rho.prior=c(0,1), NK.tuning=0.25*NK.start, d.tuning=0.25*d.start, rho.tuning=0.25*rho.start, ...)
{
	mknown <- known/N
  mK.start <- NK.start/N
  mK.tuning <- NK.tuning/N
  
  mu.values <- vector(length=size)
	sigma.values <- vector(length=size)
	mK.values <- matrix(nrow=length(indices.k), ncol=size)
	d.values <- matrix(nrow=dim(dat)[1], ncol=size)
	rho.values <- matrix(nrow=dim(dat)[2], ncol=size)
	n <- dim(dat)[1]
	
	mK.curr <- mK.start
	rho.curr <- rho.start
	d.curr <- d.start
	mu.curr <- mu.start
	sigma.curr <- sigma.start
	
	keep <- round(seq(from=burnin+1, to = burnin+iterations, length=size))
	keep.index <- 1
  
	for(t in 1:(burnin+iterations))
	{
		#cat("I am now on iteration ", t, "\n")
		
	  if(t %% 1000 == 0)
	  {
	    cat("Now on iteration ", t, "\n")
	  }
    
		d <- d.curr
		mK <- mK.curr
		rho <- rho.curr
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
				ratio <- .loglik.di.bar(dat=dat[i,], mu=mu, sigma=sigma, di=di.new, mknown=mknown, mK=mK, rho=rho) - .loglik.di.bar(dat=dat[i,], mu=mu, sigma=sigma, di=di.old, mknown=mknown, mK=mK, rho=rho)
				
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
		
		#cat("I computed d, which has a mean of ", mean(d), "\n")
		
		###########
		## update mK
		###########
		
		for(j in 1:length(indices.k))
		{
			index.this.k <- indices.k[j]
			
			mK.old <- mK[j]
			mK.new <- rnorm(1, mean=mK.old, sd=mK.tuning[j])
			
			#symmetric bounce back proposal
			while(mK.new >= 1 || mK.new <= max(dat[, index.this.k])/N)
			{
				distance.moved <- mK.new - mK.old
				if(distance.moved > 0)
				{
					max.distance.to.move <- 1 - mK.old
					mK.new <- 1 - (abs(distance.moved) - max.distance.to.move)
				} else
				{
					max.distance.to.move <- mK.old - max(dat[, index.this.k])/N
					mK.new <- max(dat[, index.this.k])/N + (abs(distance.moved) - max.distance.to.move)
				}
			} 
			
			#cat("made it out of the mk while loop \n")
			
			ratio <- .loglik.mK.bar(dat, index.this.k, n, d, mK.new, rho[index.this.k]) - .loglik.mK.bar(dat, index.this.k, n, d, mK.old, rho[index.this.k])
				
			if(is.na(ratio))
			{
				cat("I broke on mK where my K was ", index.this.k, " and my new.mK was ", mK.new, "\n")
			}
			
			accept <- min(0, ratio)
			
			logu <- log(runif(1, 0, 1))
	  
	  		if(logu < accept)
	 		{
	   			mK.curr[j] <- mK.new
	  		} else
	 		{
	    		mK.curr[j] <- mK.old
	  		}
		}
	
		
		mK <- mK.curr
		
		#cat("I computed mk, which has a mean of ", mean(mK), "\n")
		
		
		###########
		## update rho
		###########
		
		for(j in 1:dim(dat)[2])
		{
			
			rhok.old <- rho[j]
			rhok.new <- rnorm(1, mean=rhok.old, sd=rho.tuning[j])
			
			#symmetric bounce back proposal
			while(rhok.new >= rho.prior[2] || rhok.new <= rho.prior[1])
			{
				distance.moved <- rhok.new - rhok.old
				if(distance.moved > 0)
				{
					max.distance.to.move <- rho.prior[2] - rhok.old
					rhok.new <- rho.prior[2] - (distance.moved - max.distance.to.move)
				} else
				{
					max.distance.to.move <- rhok.old - rho.prior[1]
					rhok.new <- rho.prior[1] + (abs(distance.moved) - max.distance.to.move)
				}
			} 
			#cat("made it out of the rho while loop \n")
			
			mk.all <- c(mknown, mK)
	
			
			ratio <- .loglik.rhok.bar(dat, j, n, d, mk.all[j], rhok.new) - .loglik.rhok.bar(dat, j, n, d, mk.all[j], rhok.old)
				
			if(is.na(ratio))
			{
				cat("I broke on rhok where my k was ", j, " and my new.rhok was ", rhok.new, "\n")
			}
			
			accept <- min(0, ratio)
			
			logu <- log(runif(1, 0, 1))
	  
	  		if(logu < accept)
	 		{
 				rho.curr[j] <- rhok.new
	  		} else
	 		{
	    		rho.curr[j] <- rhok.old
	  		}
		}
	
		
		rho <- rho.curr
		
		#cat("I computed rho, which has a mean of ", mean(rho), "\n")
		
		###########
		## update mu
		###########

		mu.new <- rnorm(1, mean=sum(log(d))/n, sd=sigma/sqrt(n))
		
		# VERIFY  !!!!!!
		while(mu.new < mu.prior[1] || mu.new > mu.prior[2])
		{
			mu.new <- rnorm(1, mean=sum(log(d))/n, sd=sigma/sqrt(n))
		} 
		#cat("made it out of the mu while loop \n")
		
		mu.curr <- mu.new
	
		mu <- mu.curr
			
		#cat("I computed mu, which is ", mu, "\n")			
			
		###########
		## update sigma
		###########
		
		sigma.new <- sqrt(rinvgamma(1, shape=((n-1)/2), scale = 1/2*sum((log(d)-mu)^2)))
		
		while(sigma.new < sigma.prior[1] || sigma.new > sigma.prior[2])
		{
			sigma.new <- sqrt(rinvgamma(1, shape=((n-1)/2), scale = 1/2*sum((log(d)-mu)^2)))
			#cat("sigma new is ", sigma.new, "\n")
		}
		#cat("made it out of the sigma while loop \n")
		
		sigma.curr <- sigma.new
		
		#cat("I computed sigma, which is ", sigma, "\n")			
    
		if(t == keep[keep.index])
		{
		  mK.values[,keep.index] <- mK.curr
		  rho.values[,keep.index] <- rho.curr
		  d.values[,keep.index] <- d.curr
		  mu.values[keep.index] <- mu.curr
      sigma.values[keep.index] <- sigma.curr
		  keep.index <- keep.index+1   
		}
	}
	
  NK.values <- N*mK.values
	return(list(NK.values = NK.values, d.values = d.values, mu.values=mu.values, sigma.values=sigma.values, rho.values = rho.values, iterations = iterations, burnin = burnin))
}


