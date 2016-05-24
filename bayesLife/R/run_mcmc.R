if(getRversion() >= "2.15.1") utils::globalVariables("loess_sd")
data(loess_sd, envir=environment())

run.e0.mcmc <- function(sex=c("Female", "Male"), nr.chains=3, iter=160000, 
							output.dir=file.path(getwd(), 'bayesLife.output'), 
                         thin=10, replace.output=FALSE,
                         start.year=1873, present.year=2015, wpp.year=2015,
                         my.e0.file = NULL, my.locations.file = NULL, buffer.size=100,
                         a=c(13.215, 41.070, 9.235, 17.605, 2.84, 0.385),
                         #a=c(15.7669391,40.9658241,0.2107961,19.8188061,2.9306625,0.400688628),
						 delta=c(3.844, 4.035, 11.538, 5.639, 0.901, 0.4),
						 #delta=c(1.887, 1.982, 1.99, 1.949, 0.995, 0.4), 
						 tau=c(15.5976503,23.6500060,14.5056919,14.7185980,3.4514285,0.5667531), 
						 Triangle.ini = list(NULL, NULL, NULL, NULL), k.ini=NULL, z.ini=NULL,
						 lambda.ini=list(NULL, NULL, NULL, NULL), 
						 lambda.k.ini = NULL, lambda.z.ini=NULL, omega.ini = NULL,
						 Triangle.ini.low=c(10, 30, 0.1, 10), Triangle.ini.up=c(30, 50, 10, 30),
						 k.ini.low=3, k.ini.up=5, z.ini.low=0.0001, z.ini.up=0.653,
						 lambda.ini.low=c(0.01, 0.01, 0.01, 0.01), 
						 lambda.ini.up=c(0.1, 0.1, 0.1, 0.1), 
						 lambda.k.ini.low = 0.3, lambda.k.ini.up = 1, 
						 lambda.z.ini.low=1, lambda.z.ini.up=40,
						 omega.ini.low = 0.1, omega.ini.up=5, 
						 Triangle.prior.low=c(0, 0, -20, 0), Triangle.prior.up=c(100, 100, 100, 100),
						 k.prior.low=0, k.prior.up=10, z.prior.low=0, z.prior.up=0.653,
						 Triangle.c.ini.norm = list(round(Triangle.ini.low + (Triangle.ini.up - Triangle.ini.low)/2),c(2,2,2,2)), 
						 k.c.ini.norm=c(round(k.ini.low + (k.ini.up - k.ini.low)/2), 2), 
						 z.c.ini.norm=c(round(z.ini.low + (z.ini.up - z.ini.low)/2, 2), 0.2),
						 Triangle.c.prior.low=c(0, 0, -20, 0), Triangle.c.prior.up=c(100, 100, 100, 100),
						 k.c.prior.low=0, k.c.prior.up=10, z.c.prior.low=0, z.c.prior.up=0.653,
						 country.overwrites = NULL,
						 nu=4, dl.p1=9, dl.p2=9, sumTriangle.lim = c(30, 110), constant.variance=FALSE,
                         seed = NULL, parallel=FALSE, nr.nodes=nr.chains, compression.type='None',
                         auto.conf = list(max.loops=5, iter=160000, iter.incr=20000, nr.chains=3, thin=225, burnin=10000),
						 verbose=FALSE, verbose.iter = 100, ...) {
						 	
	get.init.values.between.low.and.up <- function(low, up)
		ifelse(rep(nr.chains==1, nr.chains), (low+up)/2, #seq(low, to=up, length=nr.chains)
			runif(nr.chains, low, up)
		)
		
	if(file.exists(output.dir)) {
		if(length(list.files(output.dir)) > 0 & !replace.output)
                        stop('Non-empty directory ', output.dir, 
                        ' already exists.\nSet replace.output=TRUE if you want to overwrite existing results.')
        unlink(output.dir, recursive=TRUE)
	}
    dir.create(output.dir)
    
    default.auto.conf <- formals(run.e0.mcmc)$auto.conf
	for (par in names(default.auto.conf))
		if(is.null(auto.conf[[par]])) auto.conf[[par]] <- default.auto.conf[[par]]
	auto.run <- FALSE
	if(iter == 'auto') { # defaults for auto-run (includes convergence diagnostics)
		iter <- auto.conf$iter
		nr.chains <- auto.conf$nr.chains
		auto.run <- TRUE		
	}      
	if (verbose) {
		cat('\nStarting Bayesian Hierarchical Model for Life Expectancy.\n')
		cat('=========================================================\n')
		cat('Initialize simulation -', nr.chains, 'chain(s) in total.\n')
	}
	if(!is.null(seed)) set.seed(seed)
	
	# starting values
	#==============================
	for (i in 1:4) {
		if(is.null(Triangle.ini[[i]])) 
			Triangle.ini[[i]] <- get.init.values.between.low.and.up(Triangle.ini.low[i], Triangle.ini.up[i])
		if(is.null(lambda.ini[[i]]))
			lambda.ini[[i]] <- get.init.values.between.low.and.up(lambda.ini.low[i], lambda.ini.up[i])
	}
	if(is.null(k.ini))
		k.ini <- get.init.values.between.low.and.up(k.ini.low, k.ini.up)
	if(is.null(z.ini)) 
		z.ini <- get.init.values.between.low.and.up(z.ini.low, z.ini.up)
	if(is.null(lambda.k.ini)) 
		lambda.k.ini <- get.init.values.between.low.and.up(lambda.k.ini.low, lambda.k.ini.up)
	if(is.null(lambda.z.ini)) 
		lambda.z.ini <- get.init.values.between.low.and.up(lambda.z.ini.low, lambda.z.ini.up)
	if(is.null(omega.ini)) 
		omega.ini <- get.init.values.between.low.and.up(omega.ini.low, omega.ini.up)
	sex <- substr(match.arg(sex), 1, 1)
	bayesLife.mcmc.meta <- e0.mcmc.meta.ini(sex=sex, nr.chains=nr.chains,
                                   		start.year=start.year, present.year=present.year, 
                                        wpp.year=wpp.year, my.e0.file = my.e0.file, my.locations.file=my.locations.file,
                                        output.dir=output.dir,
                                        a=a, delta=delta, tau=tau, Triangle.ini=Triangle.ini,
                                        k.ini=k.ini, z.ini=z.ini, omega.ini=omega.ini, 
                                        lambda.ini=lambda.ini, lambda.k.ini=lambda.k.ini, lambda.z.ini=lambda.z.ini,
                                        Triangle.ini.low=Triangle.ini.low, Triangle.ini.up=Triangle.ini.up, 
                                        k.ini.low=k.ini.low, k.ini.up=k.ini.up, 
                                        z.ini.low=z.ini.low, z.ini.up=z.ini.up, 
                                        Triangle.prior.low=Triangle.prior.low, Triangle.prior.up=Triangle.prior.up, 
                                        k.prior.low=k.prior.low, k.prior.up=k.prior.up, 
                                        z.prior.low=z.prior.low, z.prior.up=z.prior.up, 
                                        lambda.ini.low=lambda.ini.low, lambda.ini.up=lambda.ini.up,
                                        lambda.k.ini.low=lambda.k.ini.low, lambda.k.ini.up=lambda.k.ini.up, 
                                        lambda.z.ini.low=lambda.z.ini.low, lambda.z.ini.up=lambda.z.ini.up, 
                                        omega.ini.low=omega.ini.low, omega.ini.up=omega.ini.up, 
                                        Triangle.c.ini.norm=Triangle.c.ini.norm,
                                        k.c.ini.norm=k.c.ini.norm, z.c.ini.norm=z.c.ini.norm,
                                        Triangle.c.prior.low=Triangle.c.prior.low, Triangle.c.prior.up=Triangle.c.prior.up, 
                                        k.c.prior.low=k.c.prior.low, k.c.prior.up=k.c.prior.up, 
                                        z.c.prior.low=z.c.prior.low, z.c.prior.up=z.c.prior.up,
                                        country.overwrites=country.overwrites, 
                                        nu=nu, dl.p1=dl.p1, dl.p2=dl.p2, sumTriangle.lim=sumTriangle.lim, 
                                        constant.variance=constant.variance,
                                        buffer.size=buffer.size, compression.type=compression.type, 
                                        auto.conf=auto.conf, verbose=verbose)
    store.bayesLife.meta.object(bayesLife.mcmc.meta, output.dir)
    
    # propagate initial values for all chains if needed
    starting.values <- list()
    for (var in c('Triangle.ini', 'lambda.ini')) {
    	if(!is.list(get(var))) assign(var, list(get(var)))
    	for(i in 1:4) {
    		if (length(get(var)[[i]]) < nr.chains) {
        		if (length(get(var)[[i]]) == 1) {
            		assign(paste(var,'[[',i,']]', sep=''), rep(get(var)[[i]], nr.chains))
            	} else {
            		warning(var, '[[',i,']]', ' has the wrong length. Either 1 or ', nr.chains, 
                                ' is allowed.\nValue set to ', get(var)[[i]][1], ' for all chains.')
                    assign(paste(var,'[[',i,']]', sep=''), rep(get(var)[[i]][1], nr.chains))
               }
            }
        }
        starting.values[[var]] <- get(var)
    }
    for (var in c('k.ini', 'z.ini', 'lambda.k.ini', 'lambda.z.ini', 'omega.ini', 'iter')) {
    	if (length(get(var)) < nr.chains) {
        	if (length(get(var)) == 1) {
            	assign(var, rep(get(var), nr.chains))
            } else {
            	#stop('')
            	warning(var, ' has the wrong length. Either 1 or ', nr.chains, 
                                ' is allowed.\nValue set to ', get(var)[1], ' for all chains.')
                assign(var, rep(get(var)[1], nr.chains))
            }
        }
        if (var != 'iter') starting.values[[var]] <- get(var)
    }

	if (parallel) { # run chains in parallel
		chain.set <- bayesTFR:::bDem.performParallel(nr.nodes, 1:nr.chains, mcmc.run.chain.e0, 
                                     initfun=init.nodes.e0, meta=bayesLife.mcmc.meta, 
                                     thin=thin, iter=iter, 
                                     starting.values=starting.values,                                     
                                     verbose=verbose, verbose.iter=verbose.iter, ...)
	} else { # run chains sequentially
		chain.set <- list()
		for (chain in 1:nr.chains) {
			chain.set[[chain]] <- mcmc.run.chain.e0(chain, bayesLife.mcmc.meta, thin=thin, 
                                                iter=iter, starting.values=starting.values, 
                                                verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	names(chain.set) <- 1:nr.chains
	mcmc.set <- structure(list(meta=bayesLife.mcmc.meta, mcmc.list=chain.set), class='bayesLife.mcmc.set')
    cat('\nResults stored in', output.dir,'\n')
    
    if(auto.run) {
		diag <- try(e0.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
						thin=auto.conf$thin, burnin=auto.conf$burnin,
						verbose=verbose))
		if(auto.conf$max.loops>1) {
			for(loop in 2:auto.conf$max.loops) {
				if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
				mcmc.set <- continue.e0.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
										  parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
				diag <- try(e0.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
							thin=auto.conf$thin, burnin=auto.conf$burnin,
							verbose=verbose))
			}
		}
	}
    if (verbose) 
		cat('\nSimulation successfully finished!!!\n')
    invisible(mcmc.set)
}


mcmc.run.chain.e0 <- function(chain.id, meta, thin=1, iter=100, starting.values=NULL, verbose=FALSE, verbose.iter=10) {
	cat('\n\nChain nr.', chain.id, '\n')
    if (verbose) 
    	cat('************\n')

	this.sv <- list()
	for(var in names(starting.values)) {
		this.sv[[var]] <- if (var == 'Triangle.ini' || var == 'lambda.ini') sapply(starting.values[[var]], function(x) x[chain.id])
							else starting.values[[var]][chain.id]
	}
	mcmc <- do.call('e0.mcmc.ini', c(list(chain.id, meta, iter=iter[chain.id]), this.sv))

    if (verbose) {
        cat('Starting values:\n')
        print(unlist(mcmc[c('Triangle', 'k', 'z', 'lambda', 'lambda.k', 'lambda.z', 'omega')]))
        cat('Store initial values into ', mcmc$output.dir, '\n')
    }                
	store.e0.mcmc(mcmc, append=FALSE, flush.buffer=TRUE, verbose=verbose)

	if (verbose) 
    	cat('Start sampling -', mcmc$iter, 'iterations in total.\n')
	mcmc <- e0.mcmc.sampling(mcmc, thin=thin, verbose=verbose, verbose.iter=verbose.iter)
    return(mcmc)
}
        
continue.e0.mcmc <- function(iter, chain.ids=NULL, output.dir=file.path(getwd(), 'bayesLife.output'), 
                             parallel=FALSE, nr.nodes=NULL, auto.conf = NULL, 
                             verbose=FALSE, verbose.iter=10, ...) {
        mcmc.set <- get.e0.mcmc(output.dir)

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
                chain.list <- bayesTFR:::bDem.performParallel(nr.nodes, chain.ids, continue.e0.chain, 
                                                initfun=init.nodes.e0, mcmc.list=mcmc.set$mcmc.list, iter=iter, 
                                                verbose=verbose, verbose.iter=verbose.iter, ...)
                for (i in 1:length(chain.ids))
                        mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
        } else { # run chains sequentially
                for (chain.id in chain.ids) {
                        mcmc.set$mcmc.list[[chain.id]] <- continue.e0.chain(chain.id, mcmc.set$mcmc.list, 
                                                        iter=iter, verbose=verbose, verbose.iter=verbose.iter)
                }
        }
        cat('\n')
        if(auto.run) {
        	diag <- try(e0.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
                                    thin=auto.conf$thin, burnin=auto.conf$burnin,
                                    verbose=verbose))
			if(auto.conf$max.loops>1) {
				for(loop in 2:auto.conf$max.loops) {
					if(!inherits(diag, "try-error") && has.mcmc.converged(diag)) break
					mcmc.set <- continue.e0.mcmc(iter=auto.conf$iter.incr, output.dir=output.dir, nr.nodes=nr.nodes,
												 parallel=parallel, verbose=verbose, verbose.iter=verbose.iter, ...)
					diag <- try(e0.diagnose(sim.dir=output.dir, keep.thin.mcmc=TRUE, 
											thin=auto.conf$thin, burnin=auto.conf$burnin, verbose=verbose))
				}
			}
		}
        invisible(mcmc.set)
}

continue.e0.chain <- function(chain.id, mcmc.list, iter, verbose=FALSE, verbose.iter=10) {
        cat('\n\nChain nr.', chain.id, '\n')
        if (verbose)
                cat('************\n')
        mcmc <- mcmc.list[[chain.id]]
        mcmc$iter <- mcmc$finished.iter + iter
        if (verbose) 
                cat('Continue sampling -', iter, 'additional iterations,', mcmc$iter, 'iterations in total.\n')
        mcmc <- e0.mcmc.sampling(mcmc, thin=mcmc$thin, start.iter=mcmc$finished.iter+1, verbose=verbose, verbose.iter=verbose.iter)
        return(mcmc)
}

run.e0.mcmc.extra <- function(sim.dir=file.path(getwd(), 'bayesLife.output'), 
								countries = NULL, my.e0.file = NULL, iter = NULL,
								thin=1, burnin=2000, country.overwrites = NULL, 
								parallel=FALSE, nr.nodes=NULL, my.locations.file = NULL,
								verbose=FALSE, verbose.iter=100, ...) {
									
	mcmc.set <- get.e0.mcmc(sim.dir)
	Eini <- e0.mcmc.meta.ini.extra(mcmc.set, countries=countries, my.e0.file=my.e0.file, 
								my.locations.file=my.locations.file, burnin=burnin, 
								country.overwrites=country.overwrites, verbose=verbose)
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
		mcmc.set$mcmc.list[[chain]] <- e0.mcmc.ini.extra(mcmc.set$mcmc.list[[chain]], countries=Eini$index,
												index.replace=Eini$index.replace)
		mcthin <- max(mcthin, mcmc.set$mcmc.list[[chain]]$thin)
	}
	mcthin <- mcmc.set$mcmc.list[[1]]$thin
	total.iter <- mcmc.set$mcmc.list[[1]]$length - bayesTFR:::get.thinned.burnin(mcmc.set$mcmc.list[[1]], burnin)
	thin <- max(thin, mcthin)
	post.idx <- if (thin > mcthin) unique(round(seq(thin, total.iter, by=thin/mcthin)))
				else 1:total.iter
	if (!is.null(mcmc.set$mcmc.list[[1]]$rng.state)) .Random.seed <- mcmc.set$mcmc.list[[1]]$rng.state
	
	if (parallel) { # run chains in parallel
		if(is.null(nr.nodes)) nr.nodes<-length(chain.ids)
		chain.list <- bayesTFR:::bDem.performParallel(nr.nodes, chain.ids, e0.mcmc.run.chain.extra, 
						initfun=init.nodes.e0, mcmc.list=mcmc.set$mcmc.list, countries=Eini$index, 
						posterior.sample=post.idx, iter=iter, burnin=burnin, verbose=verbose, verbose.iter=verbose.iter, ...)
		for (i in 1:length(chain.ids))
			mcmc.set$mcmc.list[[chain.ids[i]]] <- chain.list[[i]]
	} else { # run chains sequentially
		for (chain.id in chain.ids) {
			mcmc.set$mcmc.list[[chain.id]] <- e0.mcmc.run.chain.extra(chain.id, mcmc.set$mcmc.list, 
												countries=Eini$index, posterior.sample=post.idx, iter=iter,  
												burnin=burnin, verbose=verbose, verbose.iter=verbose.iter)
		}
	}
	store.bayesLife.meta.object(meta, meta$output.dir)
	mcmc.set$meta <- meta
	cat('\n')
	invisible(mcmc.set)
}
	
e0.mcmc.run.chain.extra <- function(chain.id, mcmc.list, countries, posterior.sample, 
												iter=NULL, burnin=2000, verbose=FALSE, verbose.iter=100) {
	cat('\n\nChain nr.', chain.id, '\n')
	if (verbose)
		cat('************\n')
	mcmc <- mcmc.list[[chain.id]]
		
	if (verbose) 
		cat('MCMC sampling for additional countries and regions.\n')

	mcmc <- e0.mcmc.sampling.extra(mcmc, mcmc.list=mcmc.list, countries=countries, 
									posterior.sample=posterior.sample, 
									iter=iter, burnin=burnin, 
									verbose=verbose, verbose.iter=verbose.iter)
	return(mcmc)
}


init.nodes.e0 <- function() {
	library(bayesLife)
}
.get.Tcindex <- function(e0.matrix,  stop.if.less.than2=TRUE, cnames=NULL) {
	Tc.index <- list()
	for (country in 1:ncol(e0.matrix)) {
		Tc.index[[country]] <- which(!is.na(e0.matrix[,country]))
    	if(stop.if.less.than2 && length(Tc.index[[country]]) < 2) stop('Problem with ', cnames[country], 
    						". At least two data points must be observed.")
    }
	return(Tc.index)
}

.do.part.e0.mcmc.meta.ini <- function(data, meta) {
	nr_countries <- ncol(data$e0.matrix)
    #T_end_c <- rep(NA, nr_countries)
    Tc.index <- .get.Tcindex(data$e0.matrix, cnames=data$regions$country_name)
	T <- nrow(data$e0.matrix)
	d.ct <- loessSD <- matrix(NA, nrow=T-1, ncol=nr_countries, 
							dimnames=list(rownames(data$e0.matrix)[1:(T-1)],
									  colnames(data$e0.matrix)))
	loessSD[,] <- 1
	for(i in 2:T) {
		nisna0 <- !is.na(data$e0.matrix[i-1,])
		nisna1 <- !is.na(data$e0.matrix[i,])
		nisna2 <- nisna1 & nisna0
		if (sum(nisna2) > 0) {
			d.ct[i-1,nisna2] <- data$e0.matrix[i,nisna2] - data$e0.matrix[i-1,nisna2]
			outliers <- nisna2 & ((d.ct[i-1,] < -5) | (d.ct[i-1,] > 10))
			d.ct[i-1,outliers] <- NA
		}
		if (sum(nisna0) > 0 && !meta$constant.variance)
			loessSD[i-1,nisna0]<- loess.lookup(data$e0.matrix[i-1,nisna0])
			#loessSD[i-1,nisna0]<-sapply(data$e0.matrix[i-1,nisna0],loess.lookup)
	}
	D.supp.ct <- loessSD.suppl <- NULL
	nr_countries.suppl <- 0
	suppl <- data$suppl.data
	if(!is.null(suppl$e0.matrix)) {
		nr_countries.suppl <- ncol(suppl$e0.matrix)
		suppl$Tc.index <- .get.Tcindex(suppl$e0.matrix, stop.if.less.than2=FALSE)
		# add first time point of the observed data to get the last increment of the supplemental data
		data.suppl <- rbind(suppl$e0.matrix, data$e0.matrix[1,suppl$index.to.all.countries])
		T <- nrow(data.suppl)
		d.suppl.ct <- loessSD.suppl <- matrix(NA, nrow=T-1, ncol=nr_countries.suppl)
		for(i in 2:T) {
			nisna0 <- !is.na(data.suppl[i-1,])
			nisna1 <- !is.na(data.suppl[i,])
			nisna2 <- nisna1 & nisna0
			if (sum(nisna2) > 0) {
				d.suppl.ct[i-1,nisna2] <- data.suppl[i,nisna2] - data.suppl[i-1,nisna2]
				outliers <- nisna2 & ((d.suppl.ct[i-1,] < -5) | (d.suppl.ct[i-1,] > 10))
				d.suppl.ct[i-1,outliers] <- NA
			}
			if (sum(nisna0) > 0)
				loessSD.suppl[i-1,nisna0]<- if(meta$constant.variance) 1 else loess.lookup(data.suppl[i-1,nisna0])
				#loessSD.suppl[i-1,nisna0]<- if(meta$constant.variance) 1 else sapply(data.suppl[i-1,nisna0],loess.lookup)
		}
		suppl$nr.countries <- nr_countries.suppl
		suppl$d.ct <- d.suppl.ct
		suppl$loessSD <- loessSD.suppl
	}
	data$nr.countries <- nr_countries
	data$Tc.index <- Tc.index
	data$d.ct <- d.ct
	data$loessSD <- loessSD
	data$suppl.data <- suppl
	bounds <- .do.country.specific.ini(nr_countries, c(data, meta))
	return(c(data, bounds))
}

.do.country.specific.ini <- function(nr_countries, meta) {
	samplpars <- list()
    for (i in 1:4) {
    	samplpars[[paste('Triangle_', i, '.c.prior.low', sep='')]] <- rep(meta$Triangle.c.prior.low[i], nr_countries)
    	samplpars[[paste('Triangle_', i, '.c.prior.up', sep='')]] <- rep(meta$Triangle.c.prior.up[i], nr_countries)
    }
    samplpars$k.c.prior.low <- rep(meta$k.c.prior.low, nr_countries)
    samplpars$k.c.prior.up <- rep(meta$k.c.prior.up, nr_countries)
    samplpars$z.c.prior.low <- rep(meta$z.c.prior.low, nr_countries)
    samplpars$z.c.prior.up <- rep(meta$z.c.prior.up, nr_countries)
    for(parname in c(paste('Triangle_', 1:4, '.c.prior.low', sep=''), 
    				paste('Triangle_', 1:4, '.c.prior.up', sep=''), 'k.c.prior.low', 'k.c.prior.up', 'z.c.prior.low',
    						'z.c.prior.up'
    						))
    	names(samplpars[[parname]]) <- meta$regions$country_code

    if(!is.null(meta$country.overwrites)) {
    	for(row in 1:nrow(meta$country.overwrites)) {
    		country <- meta$country.overwrites[row, 'country_code']
    		for (col in colnames(meta$country.overwrites)) {
    			if(col == 'country_code') next
    			if(!is.element(col, names(samplpars))) {
    				warnings(col, ' is not a valid column name in country.overwrites. Column ignored.')
    				next
    			}
    			if(!is.element(as.character(country), names(samplpars[[col]]))) {
    				warnings(country, ' is not a valid country. Row ignored.')
    				break
    			}
 
    			if(!is.na(meta$country.overwrites[row,col]))
    				samplpars[[col]][as.character(country)] <- meta$country.overwrites[row,col]
    		}
    	}	
    }
	return(list(country.bounds = samplpars))
}

e0.mcmc.meta.ini <- function(sex="F", nr.chains=1, start.year=1950, present.year=2015, 
								wpp.year=2015, my.e0.file = NULL, my.locations.file = NULL,
								output.dir=file.path(getwd(), 'bayesLife.output'),
								..., verbose=FALSE) {
	mcmc.input <- c(list(sex=sex, nr.chains=nr.chains,
						start.year=start.year, present.year=present.year, 
						wpp.year=wpp.year, my.e0.file = my.e0.file,
						output.dir=output.dir), list(...))
	if(present.year-3 > wpp.year) warning("present.year is much larger then wpp.year. Make sure WPP data for present.year are available.")					
    data <- get.wpp.e0.data (sex, start.year=start.year, present.year=present.year, 
						wpp.year=wpp.year, my.e0.file = my.e0.file, 
						my.locations.file=my.locations.file, verbose=verbose)
	part.ini <- .do.part.e0.mcmc.meta.ini(data, mcmc.input)
	return(structure(c(mcmc.input, part.ini), class='bayesLife.mcmc.meta'))
}

e0.mcmc.ini <- function(chain.id, mcmc.meta, iter=100,
                     Triangle.ini = NULL, k.ini=NULL, z.ini=NULL,
				     lambda.ini=NULL, lambda.k.ini = NULL, lambda.z.ini=NULL, omega.ini = NULL,
				     verbose=FALSE) {
                                                        
	nr_countries <- mcmc.meta$nr.countries
    if (!exists(".Random.seed")) runif(1)
	mcmc <- structure(list(
						Triangle.ini=Triangle.ini,
                        k.ini=k.ini, z.ini=z.ini, omega.ini=omega.ini,
						lambda.ini=lambda.ini, lambda.k.ini=lambda.k.ini, lambda.z.ini=lambda.z.ini,
						Triangle=Triangle.ini, k=k.ini, z=z.ini, lambda=lambda.ini,
						lambda.k=lambda.k.ini, lambda.z=lambda.z.ini, omega=omega.ini,
        				output.dir=paste('mc', chain.id, sep=''), finished.iter=1, length = 1,
        				iter=iter, id=chain.id, traces=0,
        				traces.burnin=0, rng.state = .Random.seed,
        				compression.type=mcmc.meta$compression.type,
        				meta = mcmc.meta), class='bayesLife.mcmc')
    samplpars <- mcmc.meta$country.bounds
    mcmc[['Triangle.c']] <- matrix(0, ncol=nr_countries, nrow=4)
    for (i in 1:4)		
		mcmc[['Triangle.c']][i,] <- pmin(pmax(rnorm(nr_countries, mean=mcmc.meta$Triangle.c.ini.norm[[1]][i], 
										sd=mcmc.meta$Triangle.c.ini.norm[[2]][i]), 
										samplpars[[paste('Triangle_', i, '.c.prior.low', sep='')]]), 
										samplpars[[paste('Triangle_', i, '.c.prior.up', sep='')]])
	mcmc[['k.c']] <- pmin(pmax(rnorm(nr_countries, mcmc.meta$k.c.ini.norm[1], 
							sd=mcmc.meta$k.c.ini.norm[2]), samplpars$k.c.prior.low), samplpars$k.c.prior.up)
	mcmc[['z.c']] <- pmin(pmax(rnorm(nr_countries, mcmc.meta$z.c.ini.norm[1], 
							sd=mcmc.meta$z.c.ini.norm[2]), samplpars$z.c.prior.low), samplpars$z.c.prior.up)
    return(mcmc) 
}

e0.mcmc.meta.ini.extra <- function(mcmc.set, countries=NULL, my.e0.file = NULL, my.locations.file=NULL,
									burnin = 200, country.overwrites=NULL, verbose=FALSE) {
	update.regions <- function(reg, ereg, id.replace, is.new, is.old) {
		nreg <- list()
		for (name in c('code', 'area_code', 'country_code')) {
			reg[[name]][id.replace] <- ereg[[name]][is.old]
			nreg[[name]] <- c(reg[[name]], ereg[[name]][is.new])
		}
		for (name in c('name', 'area_name', 'country_name')) {
			reg[[name]][id.replace] <- as.character(ereg[[name]])[is.old]
			nreg[[name]] <- c(as.character(reg[[name]]), 
									  as.character(ereg[[name]])[is.new])
		}
		return(nreg)
	}
	update.Tc.index <- function(Tci, eTci, id.replace, is.old) {
		nTci <- Tci
		j <- length(Tci) + 1
		old.counter <- 1
		for (i in 1:length(eTci)) {
			if (is.old[i]) {
				nTci[[id.replace[old.counter]]] <- eTci[[i]]
				old.counter <- old.counter + 1
			} else {
				nTci[[j]] <- eTci[[i]]
				j <- j+1
			}
		}
		return(nTci)
	}
	update.bounds <- function(bounds, ebounds, id.replace, is.new, is.old) {
		nbounds <- list()
		for (name in c(paste('Triangle_', 1:4, '.c.prior.low', sep=''), 
    				paste('Triangle_', 1:4, '.c.prior.up', sep=''), 
    				'k.c.prior.low', 'k.c.prior.up', 'z.c.prior.low','z.c.prior.up')) {
			bounds[[name]][id.replace] <- ebounds[[name]][is.old]
			nbounds[[name]] <- c(bounds[[name]], ebounds[[name]][is.new])
		}
		return(nbounds)
	}
	meta <- mcmc.set$meta
	#create e0 matrix only for the extra countries
	e0.with.regions <- set.e0.wpp.extra(meta, countries=countries, 
									  my.e0.file = my.e0.file, my.locations.file=my.locations.file, 
									  verbose=verbose)
	if(is.null(e0.with.regions)) return(list(meta=meta, index=c()))
	meta$country.overwrites <- country.overwrites
	part.ini <- .do.part.e0.mcmc.meta.ini(e0.with.regions, meta)
	Emeta <- part.ini
						 		
	# join the new meta with the existing one
	is.old <- e0.with.regions$is_processed
	is.new <- !e0.with.regions$is_processed
	nold <- sum(is.old)
	nr_countries.all <- meta$nr.countries + Emeta$nr.countries - nold
	if (nold > 0) {
		codes.replace <- e0.with.regions$regions$country_code[is.old]
		id.replace <- unlist(sapply(codes.replace, get.country.object, meta=meta)['index',])
	} else {id.replace <- c()}
	new.meta <- list(nr.countries=nr_countries.all)
					
	for (name in c('e0.matrix', 'e0.matrix.all', 'e0.matrix.observed', 'd.ct', 'loessSD')) {
		meta[[name]][,id.replace] <- Emeta[[name]][,is.old]
		new.meta[[name]] <- cbind(meta[[name]], Emeta[[name]][,is.new])
	}
	#for (name in c('T.end.c')) {
	#	meta[[name]][id.replace] <- Emeta[[name]][is.old]
	#	new.meta[[name]] <- c(meta[[name]], Emeta[[name]][is.new])
	#}
	new.meta[['Tc.index']] <- update.Tc.index(meta$Tc.index, Emeta$Tc.index, id.replace, is.old)
	new.meta[['regions']] <- update.regions(meta$regions, Emeta$regions, id.replace, is.new, is.old)
	if(is.null(meta$country.bounds)) { # simulation was created with previous versions of bayesLife
		meta$country.bounds <- .do.country.specific.ini(meta$nr.countries, meta)$country.bounds
	}
	new.meta[['country.bounds']] <- update.bounds(meta$country.bounds, Emeta$country.bounds, id.replace, is.new, is.old)

	if(!is.null(Emeta$suppl.data$e0.matrix)) {
		suppl.id.replace <- meta$suppl.data$index.from.all.countries[id.replace]
		suppl.id.replace <- suppl.id.replace[!is.na(suppl.id.replace)]
		suppl.is.old <- which(is.old)[which(is.element(meta$suppl.data$index.from.all.countries[id.replace], suppl.id.replace))]
		suppl.old <- Emeta$suppl.data$index.from.all.countries[suppl.is.old]
		suppl.is.new <- which(is.new & !is.na(Emeta$suppl.data$index.from.all.countries))
		suppl.new <- Emeta$suppl.data$index.from.all.countries[suppl.is.new]
		for (name in c('e0.matrix', 'd.ct', 'loessSD')) {
			meta$suppl.data[[name]][,suppl.id.replace] <- Emeta$suppl.data[[name]][,suppl.old]
			new.meta$suppl.data[[name]] <- cbind(meta$suppl.data[[name]], Emeta$suppl.data[[name]][,suppl.new])
		}
		suppl.is.old.tmp <- rep(FALSE, Emeta$suppl.data$nr.countries)
		suppl.is.old.tmp[suppl.is.old] <- TRUE
		new.meta$suppl.data$Tc.index <- update.Tc.index(meta$suppl.data$Tc.index, Emeta$suppl.data$Tc.index, 
											suppl.id.replace, suppl.is.old.tmp)
		new.meta$suppl.data$regions <- update.regions(meta$suppl.data$regions, Emeta$suppl.data$regions, 
												suppl.id.replace, suppl.new, suppl.old)
		n.new <- ncol(new.meta$suppl.data$e0.matrix) - ncol(meta$suppl.data$e0.matrix)
		new.meta$suppl.data$index.from.all.countries <- meta$suppl.data$index.from.all.countries
		new.meta$suppl.data$index.to.all.countries <- meta$suppl.data$index.to.all.countries
		new.meta$suppl.data$nr.countries <- ncol(new.meta$suppl.data$e0.matrix)
		if (n.new > 0) {
			new.meta$suppl.data$index.from.all.countries <- c(new.meta$suppl.data$index.from.all.countries, rep(NA, sum(is.new)))
			new.meta$suppl.data$index.from.all.countries[meta$nr.countries + suppl.is.new] <- seq(meta$suppl.data$nr.countries + 1, 
												length=n.new)
			new.meta$suppl.data$index.to.all.countries <- c(new.meta$suppl.data$index.to.all.countries, 
											seq(meta$nr.countries+1, new.meta$nr.countries)[suppl.is.new])
		}
	}
	index <- id.replace
	if (new.meta$nr.countries > meta$nr.countries) 
		index <- c(index, seq(meta$nr.countries+1, new.meta$nr.countries))
	for (item in names(new.meta)) {
		meta[[item]] <- new.meta[[item]]
	}

	return(list(meta=meta, index=index, index.replace=id.replace))
}

e0.mcmc.ini.extra <- function(mcmc, countries, index.replace=NULL) {
	nr.countries.extra <- length(countries)
	nreplace <- length(index.replace)
	if(nreplace > 0) {
    	for (i in 1:4)		
			mcmc$Triangle.c[i,index.replace] <- pmin(pmax(rnorm(nreplace, mean=mcmc$meta$Triangle.c.ini.norm[[1]][i], 
										sd=mcmc$meta$Triangle.c.ini.norm[[2]][i]),
										mcmc$meta$Triangle.c.prior.low[i]), mcmc$meta$Triangle.c.prior.up[i])
		mcmc$k.c[index.replace] <- pmin(pmax(rnorm(nreplace, mcmc$meta$k.c.ini.norm[1], 
							sd=mcmc$meta$k.c.ini.norm[2]), mcmc$meta$k.c.prior.low), mcmc$meta$k.c.prior.up)
		mcmc$z.c[index.replace] <- pmin(pmax(rnorm(nreplace, mcmc$meta$z.c.ini.norm[1], 
							sd=mcmc$meta$z.c.ini.norm[2]), mcmc$meta$z.c.prior.low), mcmc$meta$z.c.prior.up)
	}
	samplpars <- mcmc$meta$country.bounds
	if(nr.countries.extra > nreplace) {
		nextra <- nr.countries.extra-nreplace
		eidx <- (ncol(mcmc$Triangle.c)+1):(ncol(mcmc$Triangle.c)+nextra)
		mcmc$Triangle.c <- cbind(mcmc$Triangle.c, matrix(0, ncol=nextra, nrow=4))
		for (i in 1:4)		
			mcmc$Triangle.c[i,eidx] <- pmin(pmax(rnorm(nextra, mean=mcmc$meta$Triangle.c.ini.norm[[1]][i], 
										sd=mcmc$meta$Triangle.c.ini.norm[[2]][i]),
										samplpars[[paste('Triangle_', i, '.c.prior.low', sep='')]][eidx]), 
										samplpars[[paste('Triangle_', i, '.c.prior.up', sep='')]][eidx])
		mcmc$k.c <- c(mcmc$k.c, pmin(pmax(rnorm(nextra, mcmc$meta$k.c.ini.norm[1], 
							sd=mcmc$meta$k.c.ini.norm[2]), samplpars$k.c.prior.low[eidx]), samplpars$k.c.prior.up[eidx]))
		mcmc$z.c <- c(mcmc$z.c, pmin(pmax(rnorm(nextra, mcmc$meta$z.c.ini.norm[1], 
							sd=mcmc$meta$z.c.ini.norm[2]), samplpars$z.c.prior.low[eidx]), samplpars$z.c.prior.up[eidx]))
	}
	return(mcmc)
}