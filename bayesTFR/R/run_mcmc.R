
run.tfr.mcmc <- function(nr.chains=3, iter=62000, output.dir=file.path(getwd(), 'bayesTFR.output'), 
						thin=1, replace.output=FALSE,
						# meta parameters
						start.year=1950, present.year=2015, wpp.year=2015,
						my.tfr.file = NULL, my.locations.file = NULL, buffer.size=100,
					 	U.c.low=5.5, U.up=8.8, U.width=3,
					 	mean.eps.tau0 = -0.25, sd.eps.tau0 = 0.4, nu.tau0 = 2,                                                
        				Triangle_c4.low = 1, Triangle_c4.up = 2.5,
        				Triangle_c4.trans.width=2, Triangle4.0 = 0.3, 
        				delta4.0 = 0.8, nu4 = 2,
					 	S.low=3.5, S.up=6.5, S.width=0.5,
					 	a.low=0, a.up=0.2, a.width=0.02,
					 	b.low=a.low, b.up=a.up, b.width=0.05,
					 	sigma0.low=0.01, sigma0.up=0.6, sigma0.width=0.1,
					 	sigma0.min=0.001, 
					 	const.low=0.8, const.up=2, const.width=0.3,
					 	d.low=0.05, d.up=0.5, d.trans.width=1,
					 	chi0=-1.5, psi0=0.6, nu.psi0=2,
					 	alpha0.p=c(-1, 0.5, 1.5), delta0=1, nu.delta0=2,
					 	dl.p1=9, dl.p2=9,
						# starting values (length of 1 or nr.chains)
						S.ini=NULL, a.ini=NULL, b.ini=NULL, 
					 	sigma0.ini=NULL, Triangle_c4.ini=NULL, const.ini=NULL, gamma.ini=1, 
					 	proposal_cov_gammas = NULL, # should be a list with elements 'values' and 'country_codes'
					 	seed = NULL, parallel=FALSE, nr.nodes=nr.chains, 
					 	save.all.parameters = FALSE, compression.type='None',
					 	auto.conf = list(max.loops=5, iter=62000, iter.incr=10000, nr.chains=3, thin=80, burnin=2000),
						verbose=FALSE, verbose.iter = 10, ...) {

	if(file.exists(output.dir)) {
		if(length(list.files(output.dir)) > 0 & !replace.output)
			stop('Non-empty directory ', output.dir, 
			' already exists.\nSet replace.output=TRUE if you want to overwrite existing results.')
		unlink(output.dir, recursive=TRUE)
	}
	dir.create(output.dir)
	
	default.auto.conf <- formals(run.tfr.mcmc)$auto.conf
	for (par in names(default.auto.conf))
		if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
	auto.run <- FALSE
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		iter <- auto.conf$iter
		nr.chains <- auto.conf$nr.chains
		auto.run <- TRUE		
	}
	
	if (verbose) {
		cat('\nStarting Bayesian Hierarchical Model for TFR - Phase II.\n')
		cat('========================================================\n')
		cat('Initialize simulation -', nr.chains, 'chain(s) in total.\n')
	}
	if(!is.null(seed)) set.seed(seed)
	
	# starting values (length of 1 or nr.chains)
	if (missing(S.ini) || is.null(S.ini)) 
		S.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(S.low+S.up)/2, seq(S.low, to=S.up, length=nr.chains))
	if (missing(a.ini) || is.null(a.ini)) 
		a.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(a.low+a.up)/2, seq(a.low, to=a.up, length=nr.chains))
	if (missing(b.ini) || is.null(b.ini)) 
		b.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(b.low+b.up)/2, seq(b.low, to=b.up, length=nr.chains))
	if (missing(sigma0.ini) || is.null(sigma0.ini)) 
		sigma0.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(sigma0.low+sigma0.up)/2, 
					 		seq(sigma0.low, to=sigma0.up, length=nr.chains))
	if (missing(Triangle_c4.ini) || is.null(Triangle_c4.ini)) 
		Triangle_c4.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(Triangle_c4.low+Triangle_c4.up)/2, 
					 		seq(Triangle_c4.low+0.0001, to=Triangle_c4.up-0.0001, length=nr.chains))
	if (missing(const.ini) || is.null(const.ini))  
		const.ini <- ifelse(rep(nr.chains==1, nr.chains), 
					 		(const.low+const.up)/2, 
					 		seq(const.low, to=const.up, length=nr.chains))
					 		
	bayesTFR.mcmc.meta <- mcmc.meta.ini(
						nr.chains=nr.chains,
						start.year=start.year, present.year=present.year, 
						wpp.year=wpp.year, my.tfr.file = my.tfr.file, my.locations.file=my.locations.file,
						output.dir=output.dir, phase=2,
					 	U.c.low=U.c.low, U.up=U.up, U.width=U.width,
					 	mean.eps.tau0=mean.eps.tau0, sd.eps.tau0 = sd.eps.tau0, nu.tau0 = nu.tau0,                                            
        				Triangle4.0 = Triangle4.0,  
        				Triangle_c4.low = Triangle_c4.low , Triangle_c4.up = Triangle_c4.up,
        				Triangle_c4.trans.width=Triangle_c4.trans.width,
        				delta4.0 = delta4.0, nu4=nu4,
					 	S.low=S.low, S.up=S.up, S.width=S.width,
					 	a.low=a.low, a.up=a.up, a.width=a.width,
					 	b.low=b.low, b.up=b.up, b.width=b.width,
					 	sigma0.low=sigma0.low, sigma0.up=sigma0.up, sigma0.width=sigma0.width,
					 	sigma0.min=sigma0.min, 
					 	const.low=const.low, const.up=const.up, const.width=const.width,
					 	d.low=d.low, d.up=d.up, d.trans.width=d.trans.width,
					 	chi0=chi0, psi0=psi0, nu.psi0=nu.psi0,
					 	alpha0.p = alpha0.p, delta0=delta0, nu.delta0=nu.delta0,
					 	dl.p1=dl.p1, dl.p2=dl.p2, 
					 	proposal_cov_gammas = proposal_cov_gammas,
					 	buffer.size=buffer.size, compression.type=compression.type, 
					 	auto.conf=auto.conf, verbose=verbose)
	store.bayesTFR.meta.object(bayesTFR.mcmc.meta, output.dir)
			
	# propagate initial values for all chains if needed
	for (var in c('S.ini', 'a.ini', 'b.ini', 'sigma0.ini', 'const.ini', 'gamma.ini', 'Triangle_c4.ini', 'iter')) {
		if (length(get(var)) < nr.chains) {
			if (length(get(var)) == 1) {
				assign(var, rep(get(var), nr.chains))
				} else {
				warning(var, ' has the wrong length. Either 1 or ', nr.chains, 
				' is allowed.\nValue set to ', get(var)[1], ' for all chains.')
				assign(var, rep(get(var)[1], nr.chains))
				}
			}
		}
	if (parallel) { # run chains in parallel
		chain.set <- bDem.performParallel(nr.nodes, 1:nr.chains, mcmc.run.chain, 
						initfun=init.nodes, meta=bayesTFR.mcmc.meta, 
						thin=thin, iter=iter, S.ini=S.ini, a.ini=a.ini,
                        b.ini=b.ini, sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini,
                        gamma.ini=gamma.ini, save.all.parameters=save.all.parameters, verbose=verbose, 
                        verbose.iter=verbose.iter, ...)
	} else { # run chains sequentially
		chain.set <- list()
		for (chain in 1:nr.chains) {
			chain.set[[chain]] <- mcmc.run.chain(chain, bayesTFR.mcmc.meta, thin=thin, 
					 	iter=iter, S.ini=S.ini, a.ini=a.ini, b.ini=b.ini, 
					 	sigma0.ini=sigma0.ini, Triangle_c4.ini=Triangle_c4.ini, const.ini=const.ini, 
					 	gamma.ini=gamma.ini, save.all.parameters=save.all.parameters,
					 	verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	names(chain.set) <- 1:nr.chains
	mcmc.set <- structure(list(meta=bayesTFR.mcmc.meta, mcmc.list=chain.set), class='bayesTFR.mcmc.set')
	cat('\nResults stored in', output.dir,'\n')
	
	if(auto.run) {
		diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
		if(auto.conf$max.loops>1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.tfr.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
										  parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
				diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
							thin=auto.conf$thin, burnin=auto.conf$burnin,
							verbose=verbose))
			}
		}
	}
	if (verbose)
		cat('\nSimulation successfully finished!!!\n')
	invisible(mcmc.set)
}


mcmc.run.chain <- function(chain.id, meta, thin=1, iter=100, 
							S.ini, a.ini, b.ini, sigma0.ini, Triangle_c4.ini, const.ini, gamma.ini=1,
							save.all.parameters=FALSE,
							verbose=FALSE, verbose.iter=10) {
								
	cat('\n\nChain nr.', chain.id, '\n')
    if (verbose) {
    	cat('************\n')
    	cat('Starting values:\n')
    	sv <- c(S.ini[chain.id], a.ini[chain.id], b.ini[chain.id], sigma0.ini[chain.id], Triangle_c4.ini[chain.id],
    			const.ini[chain.id], gamma.ini[chain.id])
    	names(sv) <- c('S', 'a', 'b', 'sigma0', 'Triangle_c4', 'const', 'gamma')
    	print(sv)
    }

	mcmc <- mcmc.ini(chain.id, meta, iter=iter[chain.id],
                                     S.ini=S.ini[chain.id],
                                     a.ini=a.ini[chain.id],
                                     b.ini=b.ini[chain.id],
                                     sigma0.ini=sigma0.ini[chain.id],
                                     Triangle_c4.ini=Triangle_c4.ini[chain.id],
                                     const.ini=const.ini[chain.id],
                                     gamma.ini=gamma.ini[chain.id],
                                     save.all.parameters=save.all.parameters,
                                     verbose=verbose)
	
	if (verbose) 
		cat('Store initial values into ', mcmc$output.dir, '\n')
	store.mcmc(mcmc, append=FALSE, flush.buffer=TRUE, verbose=verbose)
	
	if (verbose) 
		cat('Start sampling -', mcmc$iter, 'iterations in total.\n')
	mcmc <- tfr.mcmc.sampling(mcmc, thin=thin, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}
	
continue.tfr.mcmc <- function(iter, chain.ids=NULL, output.dir=file.path(getwd(), 'bayesTFR.output'), 
								parallel=FALSE, nr.nodes=NULL, auto.conf = NULL, verbose=FALSE, verbose.iter=10, ...) {
	mcmc.set <- get.tfr.mcmc(output.dir)

	auto.run <- FALSE
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		default.auto.conf <- mcmc.set$meta$auto.conf
		if(is.null(auto.conf)) auto.conf <- list()
		for (par in names(default.auto.conf))
			if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
		iter <- auto.conf$iter.incr
		auto.run <- TRUE
		fiter <- sapply(mcmc.set$mcmc.list, function(x) x$finished.iter)
		if (!all(fiter== fiter[1])) stop('All chains must be of the same length if the "auto" option is used.')
	}
	if (is.null(chain.ids) || auto.run) {
		chain.ids <- names(mcmc.set$mcmc.list)
	}
	if (parallel) { # run chains in parallel
		if(is.null(nr.nodes)) nr.nodes<-length(chain.ids)
		chain.list <- bDem.performParallel(nr.nodes, chain.ids, mcmc.continue.chain, 
						initfun=init.nodes, mcmc.list=mcmc.set$mcmc.list, iter=iter, verbose=verbose, 
						verbose.iter=verbose.iter, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- mcmc.continue.chain(chain.id, mcmc.set$mcmc.list, 
												iter=iter, verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	cat('\n')
	if(auto.run) {
		diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
		if(auto.conf$max.loops>1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.tfr.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
										  parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
				diag <- try(tfr.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
							thin=auto.conf$thin, burnin=auto.conf$burnin,
							verbose=verbose))
			}
		}
	}
	invisible(mcmc.set)
}
	
mcmc.continue.chain <- function(chain.id, mcmc.list, iter, verbose=FALSE, verbose.iter=10) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
	mcmc$iter <- mcmc$finished.iter + iter
	if (verbose) 
		cat('Continue sampling -', iter, 'additional iterations,', mcmc$iter, 'iterations in total.\n')

	mcmc <- tfr.mcmc.sampling(mcmc, thin=mcmc$thin, start.iter=mcmc$finished.iter+1, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}

run.tfr.mcmc.extra <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), 
								countries = NULL, my.tfr.file = NULL, iter = NULL,
								thin=1, burnin=2000, parallel=FALSE, nr.nodes=NULL, 
								my.locations.file = NULL,
								verbose=FALSE, verbose.iter=100, ...) {
									
	mcmc.set <- get.tfr.mcmc(sim.dir)
	Eini <- mcmc.meta.ini.extra(mcmc.set, countries=countries, my.tfr.file=my.tfr.file, 
												my.locations.file=my.locations.file, burnin=burnin, verbose=verbose)
	meta <- Eini$meta
	if(length(Eini$index) <= 0) {
		cat('\nNothing to be done.\n')
		return(invisible(mcmc.set))
	}
	chain.ids <- names(mcmc.set$mcmc.list)
	mcthin <- 1
	if(verbose) cat('\n')
	for (chain in chain.ids) { # update meta in each chain
		if(verbose) cat('Updating meta in chain', chain, '\n')
		mcmc.set$mcmc.list[[chain]]$meta <- meta
		mcmc.set$mcmc.list[[chain]] <- mcmc.ini.extra(mcmc.set$mcmc.list[[chain]], countries=Eini$index,
												index.replace=Eini$index.replace)
		mcthin <- max(mcthin, mcmc.set$mcmc.list[[chain]]$thin)
	}
	if(length(Eini$index_DL) <= 0) {
		cat('\nNo DL countries or regions. Nothing to be done.\n')
		store.bayesTFR.meta.object(meta, meta$output.dir)
		mcmc.set$meta <- meta
		return(invisible(mcmc.set))
	}
	mcthin <- mcmc.set$mcmc.list[[1]]$thin
	total.iter <- mcmc.set$mcmc.list[[1]]$length - get.thinned.burnin(mcmc.set$mcmc.list[[1]], burnin)
	thin <- max(thin, mcthin)
	post.idx <- if (thin > mcthin) unique(round(seq(thin, total.iter, by=thin/mcthin)))
				else 1:total.iter
	if (!is.null(mcmc.set$mcmc.list[[1]]$rng.state)) .Random.seed <- mcmc.set$mcmc.list[[1]]$rng.state
	
	if (parallel) { # run chains in parallel
		if(is.null(nr.nodes)) nr.nodes<-length(chain.ids)
		chain.list <- bDem.performParallel(nr.nodes, chain.ids, mcmc.run.chain.extra, 
						initfun=init.nodes, mcmc.list=mcmc.set$mcmc.list, countries=Eini$index_DL, 
						posterior.sample=post.idx, iter=iter, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- mcmc.run.chain.extra(chain.id, mcmc.set$mcmc.list, 
												countries=Eini$index_DL, posterior.sample=post.idx, iter=iter,  
												burnin=burnin, verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	store.bayesTFR.meta.object(meta, meta$output.dir)
	mcmc.set$meta <- meta
	cat('\n')
	invisible(mcmc.set)
}
	
mcmc.run.chain.extra <- function(chain.id, mcmc.list, countries, posterior.sample, 
												iter=NULL, burnin=2000, verbose=FALSE, verbose.iter=100) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
		
	if (verbose) 
		cat('MCMC sampling for additional countries and regions.\n')

	mcmc <- tfr.mcmc.sampling.extra(mcmc, mcmc.list=mcmc.list, countries=countries, 
									posterior.sample=posterior.sample, 
									iter=iter, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}

init.nodes <- function() {
	library(bayesTFR)
}

set.default.cltype <- function() {
#	if(!is.element(snow::getClusterOption("type"), c("MPI", "SOCK"))) 
#		snow::setDefaultClusterOptions(type="SOCK")
}

bDem.performParallel <- function(..., cltype='SOCK') {
	set.default.cltype()
	snowFT::performParallel(..., cltype=cltype)
}
