e0.mcmc.sampling <- function(mcmc, thin=1, start.iter=2, verbose=FALSE, verbose.iter=10) {
	if (!is.null(mcmc$rng.state)) .Random.seed <- mcmc$rng.state
	niter <- mcmc$iter
	T <- nrow(mcmc$meta$e0.matrix) - 1
	C <- mcmc$meta$nr.countries
	mcmc$thin <- thin	
	delta.sq <- mcmc$meta$delta^2	
	if (start.iter > niter) return(mcmc)
	
	mcenv <- new.env() # Create an environment for the mcmc stuff in order to avoid 
						# copying of the mcmc list 
    for (item in names(mcmc)) mcenv[[item]] <- mcmc[[item]]
    Triangle.prop <- rep(0,4)
    
    set.slice.sampling.width(mcenv)
    meta <- mcenv$meta
    
    dlf <- list()
    DLdata <- get.DLdata.for.estimation(meta, 1:C)
    Tm1.sum <- 0
    for(country in 1:C) Tm1.sum <- Tm1.sum + dim(DLdata[[country]])[2]
    psi.shape <- (Tm1.sum-1)/2
	recompute.par.integral <- rep(TRUE, 6)
	wpar.integral.to.mC <- sapply(1:6, compute.par.integral.to.mC, mcenv=mcenv, meta=meta, 
								lambdas.sqrt= sqrt(c(mcenv$lambda, mcenv$lambda.k, mcenv$lambda.z)), C=C)
	for(iter in start.iter:niter) {
		if(verbose.iter > 0 && (iter %% verbose.iter == 0))
			cat('\nIteration:', iter, '--', date())
		unblock.gtk('bDem.e0mcmc')

		# Update Triangle, k, z using Metropolis-Hastings sampler
		###########################################
		sum.Trkz.c <- rowSums(mcenv$Triangle.c)
		sum.Trkz.c <- c(sum.Trkz.c, sum(mcenv$k.c), sum(mcenv$z.c))
		lambdas <- c(mcenv$lambda, mcenv$lambda.k, mcenv$lambda.z)
		lambdas.sqrt <- sqrt(lambdas)
		Tr.var <- 1./(1./delta.sq + C*lambdas)
		Tr.sd <- sqrt(Tr.var)
		Tr.mean <- (meta$a/delta.sq + sum.Trkz.c*lambdas)*Tr.var
		for (i in 1:6) {
			if(recompute.par.integral[i]) {
				wpar.integral.to.mC[i] <- compute.par.integral.to.mC(i,mcenv=mcenv, meta=meta, lambdas.sqrt, C)
				recompute.par.integral[i] <- FALSE
			}
		}
		ntries <- 1
		while(ntries <= 50) {
			for (i in 1:4) {
				# Triangle is truncated normal in [0,100]
				# Triangle.prop[i] <- rnorm.trunc(mean=Tr.mean[i], sd=Tr.sd[i], 
									# low=meta$Triangle.prior.low[i], high=meta$Triangle.prior.up[i])
				# accept.prob <- ((par.integral(Triangle.prop[i], lambdas.sqrt[i], 
								# low=meta$Triangle.c.prior.low[i], up=meta$Triangle.c.prior.up[i]))^(-C))/wpar.integral.to.mC[i]
				# if (runif(1) >= accept.prob) { # not accepted
					# Triangle.prop[i] <- mcenv$Triangle[i] # replace the proposal by the old value
                    # recompute.par.integral[i] <- FALSE
                    # #cat('\nnot accepted Delta ', i, ' iter ', iter)
                # } else recompute.par.integral[i] <- TRUE
                Triangle.prop[i] <- slice.sampling(mcenv$Triangle[i], logdensity.Triangle.k.z, meta$Triangle.width[i], 
										low=meta$Triangle.prior.low[i], up=meta$Triangle.prior.up[i],
										alpha=meta$a[i], delta=meta$delta[i], par.c=mcenv$Triangle.c[i], sd=1/lambdas.sqrt[i], 
										c.low=meta$Triangle.c.prior.low[i], c.up=meta$Triangle.c.prior.up[i])
			}
			sT <- sum(Triangle.prop)
			# discard samples for which the sum(Triangle) is outside of given interval.
			if(sT <= meta$sumTriangle.lim[2] && sT >= meta$sumTriangle.lim[1]) break
			ntries <- ntries + 1
		}		
		mcenv$Triangle <- Triangle.prop
		# k is truncated normal in [0,10]
		k.prop <- rnorm.trunc(mean=Tr.mean[5], sd=Tr.sd[5], low=meta$k.prior.low, high=meta$k.prior.up)
		accept.prob <- ((par.integral(k.prop, lambdas.sqrt[5], 
								low=meta$k.c.prior.low, up=meta$k.c.prior.up))^(-C))/wpar.integral.to.mC[5]
		if (accept.prob >= 1 || runif(1) < accept.prob) {
        	mcenv$k <- k.prop
			recompute.par.integral[5] <- TRUE
		} #else {
			#recompute.par.integral[5] <- FALSE
			#cat('\nnot accepted k iter ', iter)
		#}

		# z is truncated normal in [0,1.15]
		mcenv$z <- slice.sampling(mcenv$z, logdensity.Triangle.k.z, meta$z.width, 
								low=meta$z.prior.low, up=meta$z.prior.up,
								alpha=meta$a[6], delta=meta$delta[6], par.c=mcenv$z.c, sd=1/lambdas.sqrt[6], 
								c.low=meta$z.c.prior.low, c.up=meta$z.c.prior.up)
		# z.prop <- rnorm.trunc(mean=Tr.mean[6], sd=Tr.sd[6], low=meta$z.prior.low, high=meta$z.prior.up)
		# accept.prob <- ((par.integral(z.prop, lambdas.sqrt[6], 
								# low=meta$z.c.prior.low, up=meta$z.c.prior.up))^(-C))/wpar.integral.to.mC[6]
		# if (accept.prob >= 1 || runif(1) < accept.prob) {
        	# mcenv$z <- z.prop
			# recompute.par.integral[6] <- TRUE
		# } #else {
			# #cat('\nnot accepted z iter ', iter)
			# #recompute.par.integral[6] <- FALSE
		# #}
		sum.term.for.omega <- 0
		# Update Triangle.c, k.c and z.c using slice sampling
		###########################################
		for(country in 1:C) {
			Triangle.k.z.c.update(mcenv, country, DLdata=DLdata)
			dlf[[country]] <- g.dl6(c(mcenv$Triangle.c[,country], 
								mcenv$k.c[country], mcenv$z.c[country]), 
								DLdata[[country]]['e0',], meta$dl.p1, meta$dl.p2)
			sum.term.for.omega <- sum.term.for.omega + sum(((DLdata[[country]]['dct',]-dlf[[country]])^2)/(DLdata[[country]]['loess',])^2)
		}
		# Update omega - Gibbs sampler
		###########################################
		mcenv$omega <- 1/sqrt(rgamma.ltrunc(shape=psi.shape, rate=0.5*sum.term.for.omega, low=0.01))
		#omega.update(mcenv, dlf, DLdata=DLdata)
		
		# Update lambdas using MH-algorithm
		###########################################
		for (i in 1:6) {
			if(recompute.par.integral[i]) {
				wpar.integral.to.mC[i] <- compute.par.integral.to.mC(i,mcenv=mcenv, meta=meta, lambdas.sqrt, C)
				recompute.par.integral[i] <- FALSE
			}
		}
		laccepted <- lambdas.update(mcenv, wpar.integral.to.mC, C)
		recompute.par.integral <- recompute.par.integral | laccepted
		
		################################################################### 
        # write samples simu/thin to disk
        ##################################################################
		mcenv$finished.iter <- mcenv$finished.iter+1
        mcenv$rng.state <- .Random.seed         
        if (iter %% thin == 0){
        	mcenv$length <- mcenv$length + 1
        	flush.buffer <- FALSE
            if (iter + 1 > niter) flush.buffer <- TRUE                
            store.e0.mcmc(mcenv, append=TRUE, flush.buffer=flush.buffer, verbose=verbose)
         }
	}
	resmc <- as.list(mcenv)
	class(resmc) <- class(mcmc)
	return(resmc)
}

unblock.gtk <- function(...) bayesTFR:::unblock.gtk(...)

compute.par.integral.to.mC <- function(i, mcenv, meta, lambdas.sqrt, C) {
	if(i <= 4)
		return(par.integral(mcenv$Triangle[i], lambdas.sqrt[i], 
					low=meta$Triangle.c.prior.low[i], up=meta$Triangle.c.prior.up[i])^-C)
	if(i==5)
		return(par.integral(mcenv$k, lambdas.sqrt[5], 
						low=meta$k.c.prior.low, up=meta$k.c.prior.up)^-C)
	return(par.integral(mcenv$z, lambdas.sqrt[6], 
						low=meta$z.c.prior.low, up=meta$z.c.prior.up)^-C)
}

par.integral <- function(par, sigma.inv, low=0, up=100) { 
	return(pnorm((up-par)*sigma.inv) - pnorm((low-par)*sigma.inv))
}


e0.mcmc.sampling.extra <- function(mcmc, mcmc.list, countries, posterior.sample,
											 iter=NULL, burnin=2000, 
											 verbose=FALSE, verbose.iter=100) {
	#run mcmc sampling for countries given by the index 'countries'
	niter <- mcmc$iter
	if (is.null(iter))
    	niter <- mcmc$length	
	# get values of the hyperparameters (sample from the posterior)
    hyperparameter.names <- e0.parameter.names()
    hyperparameters <- list()
    sampled.index <- sample(posterior.sample, niter, replace=TRUE)
    th.burnin <- bayesTFR:::get.thinned.burnin(mcmc, burnin)
    for (par in hyperparameter.names) {
    	hyperparameters[[par]] <- c()
    	for(mc in mcmc.list) {
    		if (bayesTFR:::no.traces.loaded(mc)  || th.burnin < mc$traces.burnin) {
    			traces <- bdem.parameter.traces(mc, par, burnin=th.burnin)
        	} else {
          		traces <- bayesTFR:::get.burned.tfr.traces(mc, par, th.burnin)
       		}
       		hyperparameters[[par]] <- rbind(hyperparameters[[par]], traces)
       	}
    	hyperparameters[[par]] <- hyperparameters[[par]][sampled.index,]
    }
	mcmc.orig <- mcmc
	mcenv <- new.env() # Create an environment for the mcmc stuff in order to avoid 
						# copying of the mcmc list 
	for (item in names(mcmc)) mcenv[[item]] <- mcmc[[item]]
	updated.var.names <- c('Triangle.c', 'k.c', 'z.c')
	DLdata <- get.DLdata.for.estimation(mcenv$meta, countries)
	set.slice.sampling.width(mcenv)
	
	for(iter in 1:niter) {
		if(verbose.iter > 0 && (iter %% verbose.iter == 0))
			cat('\nIteration:', iter, '--', date())
		unblock.gtk('bDem.e0mcmcExtra')
		# set hyperparameters for this iteration
        for (par in hyperparameter.names) {
        	if(is.null(dim(hyperparameters[[par]]))) {
        		mcenv[[par]] <- hyperparameters[[par]][iter]
        	} else {
        		mcenv[[par]] <- hyperparameters[[par]][iter,]
        	}
        }
						
		# Update Triangle.c, k.c and z.c using slice sampling
		###########################################
		for(country in countries) {
			Triangle.k.z.c.update(mcenv, country, DLdata=DLdata)
		}

		################################################################### 
        # write samples simu/thin to disk
        ##################################################################
         #update the original mcmc with the new values
         for(var in updated.var.names) {
         	mcmc.orig[[var]] <- mcenv[[var]]
         }
		 flush.buffer <- FALSE
         append <- TRUE
		 if (iter == 1) {
			append <- FALSE
			flush.buffer <- TRUE
		 } else {
			if (iter + 1 > niter) flush.buffer <- TRUE
		 }           
         store.e0.mcmc(mcmc.orig, append=append, flush.buffer=flush.buffer, countries=countries, 
         				verbose=verbose)
	}
	return(mcmc.orig)
}


get.DLdata.for.estimation <- function(meta, countries) {
	DLdata <- list()
    T.suppl.end <- if(!is.null(meta$suppl.data$e0.matrix)) nrow(meta$suppl.data$e0.matrix) else 0
    for(country in countries) {
    	idx <- which(!is.na(meta$d.ct[, country]))
    	DLdata[[country]] <- matrix(NA, nrow=3, ncol=length(idx), 
    							dimnames=list(c('e0', 'dct', 'loess'), 
    										rownames(meta$d.ct[idx,country])))
    	DLdata[[country]][1,] <- meta$e0.matrix[idx,country]
    	DLdata[[country]][2,] <- meta$d.ct[idx,country]
    	DLdata[[country]][3,] <- meta$loessSD[idx,country]
    }
    if(T.suppl.end > 0) {
    	for(country in 1:ncol(meta$suppl.data$e0.matrix)) {
    		cidx <- meta$suppl.data$index.to.all.countries[country]
    		if (!is.element(cidx, countries)) next
    		idx <- which(!is.na(meta$suppl.data$d.ct[, country]))
    		if(length(idx) <= 0) next
    		start.col <- ncol(DLdata[[cidx]]) + 1
    		DLdata[[cidx]] <- cbind(DLdata[[cidx]], 
    								matrix(NA, nrow=3, ncol=length(idx),
    										dimnames=list(rownames(DLdata[[cidx]]), 
    											rownames(meta$suppl.data$d.ct[idx,country]))))
    		end.col <- ncol(DLdata[[cidx]])
    		DLdata[[cidx]][1,start.col:end.col] <- meta$suppl.data$e0.matrix[idx,country]
       		DLdata[[cidx]][2,start.col:end.col] <- meta$suppl.data$d.ct[idx,country]
          	DLdata[[cidx]][3,start.col:end.col] <- meta$suppl.data$loessSD[idx,country]  
  		}
    }
	return(DLdata=DLdata)	
}

slice.sampling <- function(x0, fun, width,  ..., low, up, maxit=50) {
	# Slightly modified version of 
	# http://www.cs.toronto.edu/~radford/ftp/slice-R-prog (Radford M. Neal, 17 March 2008)
	gx0 <- fun(x0, ..., low=low, up=up)
	z <- gx0 - rexp(1) # slice S={x: z < gx0}
	L <- x0 - runif(1, 0, width)
	R <- L + width # should guarantee that x0 is in [L,R], even with roundoff
	#print(c(L,R,z))
	# Expand the interval until its ends are outside the slice, or until
	# the limit on steps is reached.
	J <- floor(runif(1,0,maxit))
    K <- (maxit-1) - J
    while (J>0 && L > low && fun(L,  ..., low=low, up=up)>z) {
      L <- L - width
      J <- J - 1
    }
    while (K>0 && R < up && fun(R,  ..., low=low, up=up)>z) {
      R <- R + width
      K <- K - 1
    }
    #print(c(maxit - K - J, z, L, R))
	# Shrink interval to lower and upper bounds.
	if (L<low) L <- low
  	if (R>up) R <- up
 	#if(debug) print(c('Slice sampling begin:', L, R, z, x0))
	# Sample from the interval, shrinking it on each rejection.
	i<-1
	while(i<=maxit) {
		x1 <- runif(1,L,R)
		if(z <= fun(x1,  ..., low=low, up=up)) {
			#if(debug) print(c('Slice sampling end:', L, R, x1))
			return(x1)
		}
		if (x1 < x0) L <- x1
		else R <- x1
		i <- i+1
	}
	stop('Problem in slice sampling')
}

Triangle.k.z.c.update <- function(mcmc, country, DLdata) {
	# Update Triangle pars using slice sampling
	sigmas <- 1/sqrt(mcmc$lambda)
	Triangle.c.low <- c(mcmc$meta$country.bounds$Triangle_1.c.prior.low[country], 
					  mcmc$meta$country.bounds$Triangle_2.c.prior.low[country],
					  mcmc$meta$country.bounds$Triangle_3.c.prior.low[country],
					  mcmc$meta$country.bounds$Triangle_4.c.prior.low[country])
	Triangle.c.up <- c(mcmc$meta$country.bounds$Triangle_1.c.prior.up[country], 
					  mcmc$meta$country.bounds$Triangle_2.c.prior.up[country],
					  mcmc$meta$country.bounds$Triangle_3.c.prior.up[country],
					  mcmc$meta$country.bounds$Triangle_4.c.prior.up[country])
	Triangle.prop <- rep(0,4)
	dlx <- c(mcmc$Triangle.c[,country], mcmc$k.c[country], mcmc$z.c[country])
	ntries <- 1
	while(ntries <= 50) {
		for (i in 1:4) {
			#print(c('Delta', i))
			Triangle.prop[i] <- slice.sampling(mcmc$Triangle.c[i, country],
										logdensity.Triangle.k.z.c, mcmc$meta$Triangle.c.width[i], 
										mean=mcmc$Triangle[i], 
										sd=sigmas[i], dlx=dlx,
										low=Triangle.c.low[i], 
										up=Triangle.c.up[i], par.idx=i, 
										p1=mcmc$meta$dl.p1, p2=mcmc$meta$dl.p2, omega=mcmc$omega,
										DLdata=DLdata[[country]])
			dlx[i] <- Triangle.prop[i]
		}
		sT <- sum(Triangle.prop)
		if(sT <= mcmc$meta$sumTriangle.lim[2] && sT >= mcmc$meta$sumTriangle.lim[1]) break
		dlx <- c(mcmc$Triangle.c[,country], mcmc$k.c[country], mcmc$z.c[country])
		ntries <- ntries + 1
	}
	mcmc$Triangle.c[, country] <- Triangle.prop
	#print('k')
	mcmc$k.c[country] <- slice.sampling(mcmc$k.c[country],
										logdensity.Triangle.k.z.c, mcmc$meta$k.c.width, mean=mcmc$k, 
										sd=1/sqrt(mcmc$lambda.k), dlx=dlx,
										low=mcmc$meta$country.bounds$k.c.prior.low[country], 
										up=mcmc$meta$country.bounds$k.c.prior.up[country], 
										par.idx=5, p1=mcmc$meta$dl.p1, p2=mcmc$meta$dl.p2, 
										omega=mcmc$omega, DLdata=DLdata[[country]])
	dlx[5] <- mcmc$k.c[country]
	#print('z')
	mcmc$z.c[country] <- slice.sampling(mcmc$z.c[country],
										logdensity.Triangle.k.z.c, mcmc$meta$z.c.width, mean=mcmc$z, 
										sd=1/sqrt(mcmc$lambda.z), dlx=dlx,
										low=mcmc$meta$country.bounds$z.c.prior.low[country], 
										up=mcmc$meta$country.bounds$z.c.prior.up[country], 
										par.idx=6, p1=mcmc$meta$dl.p1, p2=mcmc$meta$dl.p2, 
										omega=mcmc$omega, DLdata=DLdata[[country]])
	return()
}

logdensity.Triangle.k.z.c <- function(x, mean, sd, dlx, low, up, par.idx, p1, p2, omega, DLdata) {
	logdens <- 0.0
	res <- .C("dologdensityTrianglekz", as.double(x), as.double(mean), as.double(sd), as.double(low), 
				as.double(up), as.integer(par.idx), as.double(dlx), 
				as.double(p1), as.double(p2), as.double(DLdata['e0',]), as.integer(ncol(DLdata)), 
				as.double(DLdata['dct',]), as.double(omega*DLdata['loess',]), logdens=logdens)
	return(res$logdens)	
}

lambdas.update <- function(mcmc, wpar.integral.to.mC, C) {
	# Update lambdas using MH-algorithm
	tau.sq <- mcmc$meta$tau^2
	accepted <- rep(FALSE, 6)
	for(i in 1:4) {
		# prop <- proposal.lambda(mcmc$meta$nu, tau.sq[i], mcmc$Triangle[i], mcmc$Triangle.c[i,], 
								# mcmc$meta$nr.countries)
		# prob_accept <- ((par.integral(mcmc$Triangle[i], sqrt(prop), 
							# low=mcmc$meta$Triangle.c.prior.low[i], up=mcmc$meta$Triangle.c.prior.up[i]))^(-C))/wpar.integral.to.mC[i]
		# if (prob_accept >= 1 || runif(1) < prob_accept){
			# mcmc$lambda[i] <- prop
			# accepted[i] <- TRUE
		# }
		mcmc$lambda[i] <- slice.sampling(mcmc$lambda[i], logdensity.lambda, mcmc$meta$lambda.width[i], 
								mcmc=mcmc, low=0, up=Inf, c.low=mcmc$meta$Triangle.c.prior.low[i], c.up=mcmc$meta$Triangle.c.prior.up[i],
								Triangle=mcmc$Triangle[i], Triangle.c=mcmc$Triangle.c[i,], tau.sq=tau.sq[i])
	}
	# update lambda.k 
	prop <- proposal.lambda(mcmc$meta$nu, tau.sq[5], mcmc$k, mcmc$k.c, mcmc$meta$nr.countries)
	#lpx0 <- logdensity.lambda(mcmc$lambda.k, mcmc, mcmc$meta$k.prior.low, mcmc$meta$k.prior.up, 
	#					mcmc$k, mcmc$k.c, tau.sq[5])
	#lpx1 <- logdensity.lambda(prop, mcmc, mcmc$meta$k.prior.low, mcmc$meta$k.prior.up, mcmc$k, mcmc$k.c, tau.sq[5])
	#prob_accept <- min(exp(lpx1 - lpx0), 1)
	prob_accept <- ((par.integral(mcmc$k, sqrt(prop), 
						low=mcmc$meta$k.c.prior.low, up=mcmc$meta$k.c.prior.up))^(-C))/wpar.integral.to.mC[5]
	if (prob_accept >= 1 || runif(1) < prob_accept) {
		mcmc$lambda.k <- prop
		accepted[5] <- TRUE
	}
	# update lambda.z
	mcmc$lambda.z <- slice.sampling(mcmc$lambda.z, logdensity.lambda, mcmc$meta$lambda.z.width, 
								mcmc=mcmc, low=0, up=Inf, c.low=mcmc$meta$z.c.prior.low, c.up=mcmc$meta$z.c.prior.up,
								Triangle=mcmc$z, Triangle.c=mcmc$z.c, tau.sq=tau.sq[6])
	#prop <- proposal.lambda(mcmc$meta$nu, tau.sq[6], mcmc$z, mcmc$z.c, mcmc$meta$nr.countries)
	#lpx0 <- logdensity.lambda(mcmc$lambda.z, mcmc, mcmc$meta$z.prior.low, mcmc$meta$z.prior.up, mcmc$z, mcmc$z.c, tau.sq[6])
	#lpx1 <- logdensity.lambda(prop, mcmc, mcmc$meta$z.prior.low, mcmc$meta$z.prior.up, mcmc$z, mcmc$z.c, tau.sq[6])
	#prob_accept <- exp(lpx1 - lpx0)

	#prob_accept <- ((par.integral(mcmc$z, sqrt(prop), 
	#					low=mcmc$meta$z.c.prior.low, up=mcmc$meta$z.c.prior.up))^(-C))/wpar.integral.to.mC[6]
	#if (prob_accept >= 1 || runif(1) < prob_accept)	{
	#	mcmc$lambda.z <- prop
	#	accepted[6] <- TRUE
	#}
	return(accepted)
}

logdensity.Triangle.k.z <- function(par, low, up, alpha, delta, par.c, sd, c.low, c.up) {
	return(log(max(dnorm.trunc(par, alpha, delta, low, up), .Machine$double.xmin)) + sum(log(pmax(dnorm.trunc(par.c, par, sd, c.low, c.up),
				.Machine$double.xmin))))
}

logdensity.lambda <- function(lambda, mcmc, c.low, c.up, Triangle, Triangle.c, tau.sq, ...) {
	nu <- mcmc$meta$nu
	return(log(dgamma(lambda, nu/2, rate=tau.sq)) + 
			sum(log(pmax(dnorm.trunc(Triangle.c, mean=Triangle, sd=1/sqrt(lambda), c.low, c.up), 
				.Machine$double.xmin))))
}

proposal.lambda <- function(nu, tau.sq, Triangle, Triangle.c, C) {
	return(rgamma(1, (nu+C)/2, rate=tau.sq+0.5*sum((Triangle.c-Triangle)^2)))
}

omega.update <- function(mcmc, dlf, DLdata) {
	mcmc$omega <- slice.sampling(mcmc$omega, logdensity.omega, mcmc$meta$omega.width, 
								mcmc=mcmc, dlf=dlf, DLdata=DLdata, low=0, up=10)
	return()
}

logdensity.omega <- function(x, mcmc, dlf, DLdata, low, up) {
	log.d.cDens <- 0
	for(country in 1:mcmc$meta$nr.countries) 
		log.d.cDens <- log.d.cDens + sum(log(pmax(dnorm(DLdata[[country]]['dct',], dlf[[country]], 
										sd=x*DLdata[[country]]['loess',]),.Machine$double.xmin)))
	return(log(dunif(x, low, up)) + log.d.cDens)
}

set.slice.sampling.width <- function(mcenv){
    mcenv$meta$Triangle.c.width <- c(10, 10, 10, 10)
    mcenv$meta$Triangle.width <- c(10, 10, 10, 10)
    mcenv$meta$k.c.width <- 2
    mcenv$meta$z.c.width <- 1
    mcenv$meta$z.width <- 1
    mcenv$meta$lambda.width <- c(0.1,0.1,0.1,0.1)
    mcenv$meta$lambda.z.width <- 10
    mcenv$meta$omega.width <- 1
    return()
}

