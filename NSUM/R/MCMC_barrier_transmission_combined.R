#known includes subpopulation that will be treated as unknown in model
#matrix of tau priors give a, b terms for beta distribution for tau
#each row is for a different subpopulation
#assumes unknown subpopulations that have tau's are last in "known"
.simulate.comb <- function(n, known, unknown, N, mu=5, sigma=1, rho=rep(0.1, length(c(known, unknown))), tauK=rep(1, length(unknown)), ...)
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
	
	#draw tau_K for unknown subpopulations
	tau <-c(rep(1, length(known)), tauK)
	
	y <- matrix(nrow=n, ncol=length(known.all))
	for(i in 1:n)
	{
		di <- d[i]
		for(k in 1:length(known.all))
		{
			y[i,k] <- rbinom(1, size=as.integer(di), prob=tau[k]*proportions[i,k])
		}
	}
	
	return(list(y=y, d=d))
}

#function for mK terms of log posterior 
#mK - estimated value for m_K
#n - is number of respondents in survey
#rhoK - rho for a given mK
#index.k - for term being updated
.loglik.mK.comb <- function(dat, index.k, n, d, mK, rhoK, q)
{
	q.relevant <- q[,index.k]
	
	rhoK.term <- 1/rhoK - 1 #to be computationally faster
	mK.rhoK.term <- mK*rhoK.term
	
	sum.qik <- (mK.rhoK.term - 1)*sum(log(q.relevant))
	sum.minusqik <-(rhoK.term - mK.rhoK.term - 1)*sum(log(1-q.relevant))
	beta.term <- lbeta(mK.rhoK.term, rhoK.term - mK.rhoK.term)
	
	loglik <- sum.qik + sum.minusqik - n*beta.term - log(mK)
	
	return(loglik)
}


#function for rho_k terms of log posterior 
#rhok - estimated value for rho_k
#n - is number of respondents in survey
#mk - single m_k value (rho_k matches for the k)
#index.k - index of K being updated
.loglik.rhok.comb <- function(dat, index.k, n, d, mk, rhok, q)
{
	q.relevant <- q[,index.k]
	
	#to be computationally faster
	rhok.term <- 1/rhok - 1 
	mk.rhok.term <- mk*rhok.term
	
	sum.qik <- (mk.rhok.term - 1)*sum(log(q.relevant))
	sum.minusqik <- (rhok.term - mk.rhok.term - 1)*sum(log(1-q.relevant))
	beta.term <- lbeta(mk.rhok.term, rhok.term - mk.rhok.term)
	
	loglik <- sum.qik + sum.minusqik - n*beta.term
	
	return(loglik)	
}



#dat - all data
#d - all d
#index.i - relevant index i
#index.k - relevant index k
#mk.all - all mk (i.e. c(known/N, mK))
#tauk.all - all tauK (1's for all known tauK)
#qik - relevant qik term

.loglik.qik.comb <- function(dat, d, index.i, index.k, mk.all, tauk.all, rhok, qik)
{
	dat.relevant <- dat[index.i, index.k]
	d.relevant <- d[index.i]
	mk.relevant <- mk.all[index.k]
	rhok.relevant <- rhok[index.k]
	tauk.relevant <- tauk.all[index.k]
	
	rhok.term <- 1/rhok.relevant - 1
	mk.rhok.term <- mk.relevant * rhok.term
	
	qik.term <- (dat.relevant + mk.rhok.term - 1) * log(qik)
	qik.tauk.term <- (d.relevant - dat.relevant) * log(1-tauk.relevant*qik)
	qik.minus.term <- (rhok.term - mk.rhok.term - 1)*log(1-qik)
	
	loglik <- qik.term + qik.tauk.term + qik.minus.term
	
	return(loglik)
}

#function for tau_k terms of log posterior 
#tauK - estimated value for rho_K
#n - is number of respondents in survey
#index.k - index of K being updated
# dat - all dat
# q - all q_ik
#tauK.prior - relevant prior
.loglik.tauK.comb <- function(dat, index.k, n, d, q, tauK, tauK.prior)
{
	aK <- tauK.prior[1]
	bK <- tauK.prior[2]
	
	q.relevant <- q[,index.k]
	dat.relevant <- dat[,index.k]
	
	tauk.term <- sum(dat.relevant + aK - 1)*log(tauK)
	minus.tauk.qik.term <- sum((d - dat.relevant)*log(1-tauK*q.relevant))
	minus.tauk.term <- n*(bK-1)*log(1-tauK)
	
	
	loglik <- tauk.term + minus.tauk.qik.term + minus.tauk.term
	
	return(loglik)	
}




#function for d_i log posterior terms
#takes into dat only for individual i
#takes in pik only for individual i
#takes in tauk.all - 1's for all known subpops, tauK for K unknown
.loglik.di.comb <- function(dat, mu, sigma, di, qik, tauk.all)
{	
	sum.binom <- sum(lchoose(di, dat))
	
	tau.qik.term <- di * sum(log(1-qik*tauk.all))
	
	loglik <- -log(di) - (log(di) - mu)^2 / (2*sigma^2) + sum.binom + tau.qik.term
	
	return(loglik)
}


#mKstart, mK.tuning, tauK.start, tauK.tuning, indices.k should all be same length (length of # of N_K being sampled)
#mknown - known values divided by N
#q - n x K matrix at each iteration of probability individual i knows someone in group k
#mu.prior, sigma.prior - length 2 vector having lower and upper limits of uniform prior
#tauk.priors - matrix with row for each unknow, 1st column is a_k, 2nd column is b_k
  
.mcmc.comb <- function(dat, known, N, indices.k, iterations, burnin, size, NK.start=killworth.start(dat,known,N)$NK.start, d.start=killworth.start(dat,known,N)$d.start, mu.start=killworth.start(dat,known,N)$mu.start, sigma.start=killworth.start(dat,known,N)$sigma.start, rho.start=rep(0.1,dim(dat)[2]), tauK.start=rep(0.5,length(indices.k)), q.start=matrix(c(known,NK.start)/N,dim(dat)[1],dim(dat)[2]), mu.prior=c(3,8), sigma.prior=c(1/4,2), rho.prior=c(0,1), tauK.priors=matrix(rep(c(1,1),each=length(indices.k)),length(indices.k),2), NK.tuning=0.25*NK.start, d.tuning=0.25*d.start, rho.tuning=0.25*rho.start, tauK.tuning=0.25*tauK.start, q.tuning=0.25*q.start, ...)
{
  mknown <- known/N
  mK.start <- NK.start/N
  mK.tuning <- NK.tuning/N
  
  n <- dim(dat)[1]
  num.known <- dim(dat)[2] - length(indices.k)
  mu.values <- vector(length=size)
  sigma.values <- vector(length=size)
  mK.values <- matrix(nrow=length(indices.k), ncol=size)
  d.values <- matrix(nrow=n, ncol=size)
  rho.values <- matrix(nrow=dim(dat)[2], ncol=size)
  tauK.values <- matrix(nrow=length(indices.k), ncol=size)
  q.values <- array(dim=c(n, dim(dat)[2], size))
  
  mK.curr <- mK.start
  rho.curr <- rho.start
  tauK.curr <- tauK.start
  q.curr <- q.start
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
    tauK <- tauK.curr
    q <- q.curr
    mu <- mu.curr
    sigma <- sigma.curr
    
    #############
    # compute d
    #############
    
    
    tauk.all <- c(rep(1, num.known), tauK)
    
    for(i in 1:n)
    {
      di.old <- d[i]
      di.new <- rnorm(1, mean=di.old, sd=d.tuning[i])
      
      if(di.new < max(dat[i,]))
      {
        d.curr[i] <- di.old
      } else
      {				
        ratio <- .loglik.di.comb(dat=dat[i,], mu=mu, sigma=sigma, di=di.new, qik = q[i,], tauk.all=tauk.all) - .loglik.di.comb(dat=dat[i,], mu=mu, sigma=sigma, di=di.old, qik = q[i,], tauk.all=tauk.all)
        
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
      
      ratio <- .loglik.mK.comb(dat, index.this.k, n, d, mK.new, rho[index.this.k], q=q) - .loglik.mK.comb(dat, index.this.k, n, d, mK.old, rho[index.this.k], q=q)
      
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
    
    mk.all <- c(mknown, mK)
    
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
      
      ratio <- .loglik.rhok.comb(dat, j, n, d, mk.all[j], rhok.new, q) - .loglik.rhok.comb(dat, j, n, d, mk.all[j], rhok.old, q)
      
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
    ## update tauK
    ###########
    
    for(j in 1:length(indices.k))
    {
      index.this.k <- indices.k[j]
      
      tauK.old <- tauK[j]
      tauK.new <- rnorm(1, mean=tauK.old, sd=tauK.tuning[j])
      
      #symmetric bounce back proposal
      while(tauK.new >= 1 || tauK.new <= 0)
      {
        distance.moved <- tauK.new - tauK.old
        if(distance.moved > 0)
        {
          max.distance.to.move <- 1 - tauK.old
          tauK.new <- 1 - (abs(distance.moved) - max.distance.to.move)
        } else
        {
          max.distance.to.move <- tauK.old - 0
          tauK.new <- 0 + (abs(distance.moved) - max.distance.to.move)
        }
      } 
      
      #cat("made it out of the mk while loop \n")
      
      ratio <- .loglik.tauK.comb(dat, index.this.k, n, d, q, tauK.new, tauK.priors[j,]) - .loglik.tauK.comb(dat, index.this.k, n, d, q, tauK.old, tauK.priors[j,]) 
      
      if(is.na(ratio))
      {
        cat("I broke on tauK where my K was ", index.this.k, " and my new.tauK was ", tauK.new, "\n")
      }
      
      accept <- min(0, ratio)
      
      logu <- log(runif(1, 0, 1))
      
      if(logu < accept)
      {
        tauK.curr[j] <- tauK.new
      } else
      {
        tauK.curr[j] <- tauK.old
      }
    }
    
    
    tauK <- tauK.curr
    
    tauk.all <- c(rep(1, num.known), tauK)
    
    #cat("I computed tauk, which has a mean of ", mean(tauK), "\n")
    
    
    ###########
    ## update q
    ###########
    
    for(i in 1:n)
    {
      for(k in 1:dim(dat)[2])
      {
        qik.old <- q[i,k]
        qik.new <- rnorm(1, mean=qik.old, sd=q.tuning[i,k])
        
        #symmetric bounce back proposal
        while(qik.new >= 1 || qik.new <= 0)
        {
          distance.moved <- qik.new - qik.old
          if(distance.moved > 0)
          {
            max.distance.to.move <- 1 - qik.old
            qik.new <- 1 - (abs(distance.moved) - max.distance.to.move)
          } else
          {
            max.distance.to.move <- qik.old - 0
            qik.new <- 0 + (abs(distance.moved) - max.distance.to.move)
          }
        } 
        
        #cat("made it out of the qik while loop \n")
        
        ratio <- .loglik.qik.comb(dat, d, i, k, mk.all, tauk.all, rho, qik.new) - .loglik.qik.comb(dat, d, i, k, mk.all, tauk.all, rho, qik.old)
        
        if(is.na(ratio))
        {
          cat("I broke on qik for i ", i, " and k ", k, " where my new qik was ", qik.new, "\n")
        }
        
        accept <- min(0, ratio)
        
        logu <- log(runif(1, 0, 1))
        
        if(logu < accept)
        {
          q.curr[i,k] <- qik.new
        } else
        {
          q.curr[i,k] <- qik.old
        }
      }
    }
    
    
    q <- q.curr
    
    
    #cat("I computed q, which has a mean of ", mean(q), "\n")
    
    ###########
    ## update mu
    ###########
    
    mu.new <- rnorm(1, mean=sum(log(d))/n, sd=sigma/sqrt(n))
    
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
      mu.values[keep.index] <- mu.curr
      sigma.values[keep.index] <- sigma.curr
      mK.values[,keep.index] <- mK.curr
      d.values[,keep.index] <- d.curr
      rho.values[,keep.index] <- rho.curr
      tauK.values[,keep.index] <- tauK.curr
      q.values[,,keep.index] <- q.curr  
      keep.index <- keep.index+1   
    }
  }
  
  NK.values <- N*mK.values
  return(list(NK.values = NK.values, d.values = d.values, mu.values=mu.values, sigma.values=sigma.values, rho.values = rho.values, tauK.values = tauK.values, q.values = q.values, iterations = iterations, burnin = burnin))
}