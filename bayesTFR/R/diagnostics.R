tfr.raftery.diag <- function(mcmc=NULL, 
							 sim.dir=file.path(getwd(), 'bayesTFR.output'),
							 burnin=0, country=NULL,
							 par.names = tfr.parameter.names(trans=TRUE),
							 par.names.cs =tfr.parameter.names.cs(trans=TRUE, back.trans=FALSE),
							 country.sampling.prop=1,
							 verbose=TRUE, ...
							 ) {
	is.error <- function(rd) {
		if(is.null(dim(rd[[1]]$resmatrix))) {
			print(rd)
			return(TRUE)
		}
		return(FALSE)
	}
	if(verbose) cat('\nComputing Raftery Diagnostics ... this can take a while ...\n')
	if (is.null(mcmc)) {
		mcmc.set <- get.tfr.mcmc(sim.dir=sim.dir, low.memory=TRUE)
	} else {
		mcmc.set <- mcmc
	}
	gui.option.name <- paste('bDem', strsplit(class(mcmc.set), '.', fixed=TRUE)[[1]][1], 'diagnose', sep='.')
	gui.option.name.status <- paste(gui.option.name, 'status', sep='.')
	gui.options <- list()
	result.025 <- result.975 <- burnin.025 <- burnin.975 <- rd.025<- thin.ind.025 <- thin.ind.975 <- NULL
	if (!is.null(par.names)) {
		coda.mc <- coda.list.mcmc(mcmc.set, country=country, rm.const.pars=TRUE, low.memory=TRUE,
						par.names=par.names, par.names.cs=NULL, burnin=burnin, ...
						)
		gui.options[[gui.option.name.status]] <- 'processing country-independent pars (q=0.025)'
		unblock.gtk(gui.option.name, gui.options)
		if(verbose) cat('\t\tProcessing raftery.diag(..., r=0.0125, q=0.025) for country-independent parameters\n')
		rd.025 <- raftery.diag(coda.mc, r=0.0125, q=0.025)
		if(is.error(rd.025)) return()
		thin.ind.025 <- diag.thin.indep(coda.mc, q=0.025)
		colnames(thin.ind.025) <- rownames(rd.025[[1]]$resmatrix)
		gui.options[[gui.option.name.status]] <- 'processing country-independent pars (q=0.975)'
		unblock.gtk(gui.option.name, gui.options)
		if(verbose) cat('\t\tProcessing raftery.diag(..., r=0.0125, q=0.975) for country-independent parameters\n')
		rd.975 <- raftery.diag(coda.mc, r=0.0125, q=0.975)
		if(is.error(rd.975)) return()
		thin.ind.975 <- diag.thin.indep(coda.mc, q=0.975)
		colnames(thin.ind.975) <- rownames(rd.025[[1]]$resmatrix)
		result.025 <- result.975 <- burnin.025 <- burnin.975 <- matrix(NA, 
											nrow=length(rd.025), ncol=dim(rd.025[[1]]$resmatrix)[1],
											dimnames=list(c(), rownames(rd.025[[1]]$resmatrix)))
		for(i in 1:length(rd.025)) {
			result.025[i,] <- rd.025[[i]]$resmatrix[,'N']
			result.975[i,] <- rd.975[[i]]$resmatrix[,'N']
			burnin.025[i,] <- rd.025[[i]]$resmatrix[,'M']
			burnin.975[i,] <- rd.975[[i]]$resmatrix[,'M']
		}
	}
	c.index <- NULL
	if(!is.null(par.names.cs)) {
		if (!is.null(country)) {
			country.obj <- get.country.object(country, mcmc.set$meta)
			c.index <- country.obj$index
		} else {
			c.index <- get.countries.index(mcmc.set$meta)
			if(country.sampling.prop < 1) 
				c.index <- sort(sample(c.index, size=round(length(c.index)*country.sampling.prop,0)))
		}
		if(verbose) 
			cat('\t\tProcessing raftery.diag for country ')
		country.counter <- 0
		status.for.gui <- paste('out of', length(c.index), 'countries.')
		for(country.idx in c.index) {
			if(getOption(gui.option.name, default=FALSE)) {
				# This is to unblock the GUI, if the run is invoked from bayesDem
				# and pass info about its status
				# In such a case the gtk libraries are already loaded
				country.counter <- country.counter + 1
				gui.options[[gui.option.name.status]] <- paste('finished', country.counter, status.for.gui)
				unblock.gtk(gui.option.name, gui.options)
			}
			country.obj <- get.country.object(country.idx, mcmc.set$meta, index=TRUE)
			coda.mc.cs <- coda.list.mcmc(mcmc.set, country=country.obj$code, 
									rm.const.pars=TRUE, low.memory=TRUE,
									par.names = NULL,
									par.names.cs=par.names.cs, burnin=burnin, ...
									)
			if(verbose) 
				cat(country.idx, ', ')
			rd.025 <- raftery.diag(coda.mc.cs, r=0.0125, q=0.025)
			if(is.error(rd.025)) return()
			thin.ind.025 <- cbind(thin.ind.025, diag.thin.indep(coda.mc.cs, q=0.025))
			npar <- dim(rd.025[[1]]$resmatrix)[1]
			ncols <- ncol(thin.ind.025)
			colnames(thin.ind.025)[(ncols-npar+1):ncols] <- rownames(rd.025[[1]]$resmatrix)
			rd.975 <- raftery.diag(coda.mc.cs, r=0.0125, q=0.975)
			if(is.error(rd.975)) return()
			thin.ind.975 <- cbind(thin.ind.975, diag.thin.indep(coda.mc.cs, q=0.975))
			colnames(thin.ind.975)[(ncols-npar+1):ncols] <- rownames(rd.025[[1]]$resmatrix)
			
			m <- matrix(NA, nrow=length(rd.025), ncol=npar, 
							dimnames=list(c(), rownames(rd.025[[1]]$resmatrix)))
			result.025 <- cbind(result.025, m)
			result.975 <- cbind(result.975, m)
			burnin.025 <- cbind(burnin.025, m)
			burnin.975 <- cbind(burnin.975, m)
			for(i in 1:length(rd.025)) {
				idx <- (ncol(result.025)-npar+1):ncol(result.025)
				result.025[i, idx] <- rd.025[[i]]$resmatrix[,'N']
				result.975[i, idx] <- rd.975[[i]]$resmatrix[,'N']
				burnin.025[i, idx] <- rd.025[[i]]$resmatrix[,'M']
				burnin.975[i, idx] <- rd.975[[i]]$resmatrix[,'M']
			}
			gc()
		}
		if(verbose) 
			cat('\n')		
	}
	if(is.null(rd.025)) {
		if(verbose) cat('\t\tNothing to be done.\n')
		return()
	}
	if(verbose) cat('\t\tProcessing results ...')
	nr.chains <- length(rd.025)
	nr.par <- ncol(result.025)
	par.names.all <- colnames(result.025)
	if (is.null(par.names.all)) { # probably an error occured in raftery.diag
		print(rd.025)
		return()
		}

	colidx.cind <- get.full.par.names(par.names, par.names.all, index=TRUE)
	colidx.cdep <- get.full.par.names.cs(par.names.cs, par.names.all, index=TRUE)
	
	lcolidx.cind <- length(colidx.cind)
	lcolidx.cdep <- length(colidx.cdep)
	N.country.ind <- N.country.dep <- notconv1 <- notconv2 <- notconvinchain1 <- notconvinchain2 <- NULL
	iter <- nr.chains*(mcmc.set$mcmc.list[[1]]$finished.iter - burnin)
	for(chain in 1:nr.chains) {
		N <- rbind(result.025[chain,], result.975[chain,])
		where.larger <- N > iter
		notconv1 <- rbind(notconv1, 
								  cbind(parameter.name=par.names.all[where.larger[1,]], 
										chain.id=rep(chain, sum(where.larger[1,])),
										N=N[1,where.larger[1,]])
								  )
		notconv2 <- rbind(notconv2, 
								  cbind(parameter.name=par.names.all[where.larger[2,]], 
										chain.id=rep(chain, sum(where.larger[2,])),
										N=N[2,where.larger[2,]])
								  )
		where.larger.inchain <- N > (mcmc.set$mcmc.list[[chain]]$finished.iter - burnin)
		notconvinchain1 <- rbind(notconvinchain1, 
								  cbind(parameter.name=par.names.all[where.larger.inchain[1,]], 
										chain.id=rep(chain, sum(where.larger.inchain[1,])),
										N=N[1,where.larger.inchain[1,]])
								  )
		notconvinchain2 <- rbind(notconvinchain2, 
								  cbind(parameter.name=par.names.all[where.larger.inchain[2,]], 
										chain.id=rep(chain, sum(where.larger.inchain[2,])),
										N=N[2,where.larger.inchain[2,]])
								  )

		N.country.ind <- rbind(N.country.ind,
								cbind(parameter.name=par.names.all[colidx.cind],
									  chain.id=rep(chain, lcolidx.cind),
									  N0.025=N[1,colidx.cind],
									  N0.975=N[2,colidx.cind]))
		N.country.dep <- rbind(N.country.dep,
								cbind(parameter.name=par.names.all[colidx.cdep],
									  chain.id=rep(chain, lcolidx.cdep),
									  N0.025=N[1,colidx.cdep],
									  N0.975=N[2,colidx.cdep]))
	}
	notconv <- list(notconv1, notconv2)
	notconvinchain <- list(notconvinchain1, notconvinchain2)
	names(notconv) <- names(notconvinchain) <- c('0.025', '0.975')
	# change data types
	for (i in 1:2) {
		if (nrow(notconv[[i]]) > 0) {
			notconv[[i]] <- data.frame(notconv[[i]], row.names=1:nrow(notconv[[i]]))
		}else{notconv[[i]] <- NA}
		if (nrow(notconvinchain[[i]]) > 0) {
			notconvinchain[[i]] <- data.frame(notconvinchain[[i]], row.names=1:nrow(notconvinchain[[i]]))
		}else{notconvinchain[[i]] <- NA}
	}
	if (!is.null(par.names)) {
		N.country.ind <- data.frame(N.country.ind, row.names=1:nrow(N.country.ind))
	} else {N.country.ind <- NULL}
	if (!is.null(par.names.cs)) {
		N.country.dep <- data.frame(N.country.dep, row.names=1:nrow(N.country.dep))
	} else {N.country.dep <- NULL}
	
	#compute medians over chains
	s <- rbind(apply(result.025, 2, median), apply(result.975, 2, median))
	bi <- rbind(apply(burnin.025, 2, median), apply(burnin.975, 2, median))
	thin.ind <- rbind(apply(thin.ind.025, 2, median), apply(thin.ind.975, 2, median))

	colnames(s) <- par.names.all

	full.par.names <- get.full.par.names(par.names, par.names.all)
	full.par.names.cs <- get.full.par.names.cs(par.names.cs, par.names.all)
	if (is.null(country) && !is.null(full.par.names.cs)) {
		short.par.names.cs <- strsplit(full.par.names.cs, "_c.*") # removes the '_c*' postfix
		short.par.names.cs <- unique(short.par.names.cs)
	} else short.par.names.cs <- full.par.names.cs
	qcs <- qbics <- matrix(NA, nrow=2, ncol=length(short.par.names.cs)+length(full.par.names),
						dimnames=list(c(), c(short.par.names.cs, full.par.names)))
	if (is.null(country)) {
		# compute maximum for U, d, and gammat_1/2/3 over countries
		for (pname in short.par.names.cs) {
			colidx <- grep(paste('^', pname, '_', sep=''), par.names.all)
			for (j in 1:2) {
				qcs[j,pname] <- max(s[j,colidx])
				qbics[j,pname] <- max(bi[j,colidx])
			}
		}
	} else {
		qcs[,full.par.names.cs] <- s[, full.par.names.cs]
		qbics[,full.par.names.cs] <- bi[, full.par.names.cs]
	}
	Nmed.cs <- s[, full.par.names.cs, drop=FALSE]
	qcs[,full.par.names] <- s[, full.par.names]
	qbics[,full.par.names] <- bi[, full.par.names]
	if(verbose) cat(' done.\n')
	return(list(Nmedian=round(qcs,0),
				burnin=round(qbics,0),
				not.converged.parameters=notconv,
				not.converged.inchain.parameters=notconvinchain,
				N.country.indep=N.country.ind,
				N.country.spec=N.country.dep,
				Nmedian.country.spec=round(Nmed.cs,0),
				thin.ind=list('0.025'=thin.ind.025, '0.975'=thin.ind.975, median=thin.ind),
				nr.countries=c(used=length(c.index), total=get.nr.countries.est(mcmc.set$meta))))
}

process.not.converged.parameters <- function(diag, iter) {
	diff <-apply(diag$Nmedian, 2, max) - iter
	posdiff <- diff[diff > 0]
	not.converged <- names(posdiff)
	lnot.converged <- length(not.converged)
	N <- matrix(0, nrow=lnot.converged, ncol=2, dimnames=list(not.converged,
					c('Total iterations needed', 'Remaining iterations')))
	if(lnot.converged == 0) return(N)
	for (i in 1:2) {
		is.contained <- is.element(not.converged, 
							colnames(diag$Nmedian))
		parnames <- not.converged[is.contained]
		for (parname in parnames) {
			N[parname,1] <- max(N[parname,1], diag$Nmedian[i, parname])
		}
	}
	N[,2] <- pmax(N[,1] - iter, 0)
	return(N)
}
	
.do.diagnose <- function(type, class.name, sim.dir, thin=80, burnin=2000, express=FALSE, 
							country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
	get.country.name <- function(par) {
		cindex <- strsplit(par, '_c')[[1]]
		cindex <- as.numeric(cindex[length(cindex)])
		return(get.country.object(cindex, meta=mcmc.set$meta)$name)
	}

	mcmc.set <- do.call(paste('get.', type, '.mcmc', sep=''), 
						list(sim.dir=sim.dir, low.memory=TRUE))
	if(is.null(mcmc.set))
		stop('No valid simulation in ', sim.dir)
	iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
	if(iter <= 0) stop('0 number of iterations. Check the value of burnin.')
	if(iter/thin < 1) stop(paste('Value of thin is too high (', thin, 
							') given the total number of iterations (', iter, ').', sep=''))
	#run raftery.diag on country-independent parameters
	diag.procedure <- paste(type, '.raftery.diag', sep='')
	raftery.diag.res <- do.call(diag.procedure, 
								list(mcmc.set, par.names.cs=NULL, thin=thin,
								burnin=burnin, verbose=verbose))
	if (is.null(raftery.diag.res)) stop(paste('Problem in', diag.procedure))
	raftery.diag.res.cs <- NULL
	if(!is.null(country.sampling.prop)) express <- FALSE
	if(!express &&((!is.null(country.sampling.prop) && (country.sampling.prop>0)) 
						|| is.null(country.sampling.prop))) {
		#run raftery.diag on country-specific parameters
		raftery.diag.res.cs <- do.call(diag.procedure, list(mcmc.set, 
								par.names = NULL,
								thin=thin,  burnin=burnin,
								country.sampling.prop=if(is.null(country.sampling.prop)) 1 else country.sampling.prop,
								verbose=verbose
								))
		if (is.null(raftery.diag.res.cs)) stop(paste('Problem in', diag.procedure))
	}
	status <- c(red=FALSE, green=FALSE)
	res <- NULL
	lres.cind <- 0
	if (!is.null(raftery.diag.res)) {
		res <- process.not.converged.parameters(raftery.diag.res, iter)
		lres.cind <- nrow(res)
	}
	if (!is.null(raftery.diag.res.cs)) {
		res <- rbind(res, process.not.converged.parameters(raftery.diag.res.cs, iter))
	} 
	to.run <- 0
	if(nrow(res) > 0) {
		max.idx <- which.max(res[,2])
		to.run <- res[max.idx,2]
	}
	if(to.run <= 0) status['green'] <- TRUE
	else status['red'] <- TRUE
	nr.countries <- if(!is.null(raftery.diag.res.cs)) raftery.diag.res.cs$nr.countries 
					else c(used=0, total=get.nr.countries.est(mcmc.set$meta))
	use.nr.traj <- floor(iter/thin)
	thinned.mcmc <- NULL
	if(keep.thin.mcmc) {
		thinned.mcmc <- do.call(paste('get.thinned.', type, '.mcmc', sep=''), 
										list(mcmc.set, thin=thin, burnin=burnin))
		if(is.null(thinned.mcmc) || thinned.mcmc$meta$parent.iter < iter) 
			thinned.mcmc <- do.call(paste('create.thinned.', type, '.mcmc', sep=''), 
										list(mcmc.set, thin=thin, burnin=burnin, verbose=verbose))
	}
	diag <- structure(list(result=res,
					lresult.country.independent=lres.cind,
					country.independent=raftery.diag.res, 
					country.specific=raftery.diag.res.cs, 
					iter.needed=to.run,
					iter.total=iter, use.nr.traj=use.nr.traj,
					status=status,
					mcmc.set=mcmc.set,
					thin.mcmc=thinned.mcmc,
					burnin=burnin,
					thin=thin,
					express=express, 
					nr.countries=nr.countries),
				 class=class.name)
	if(verbose) summary(diag)
	save.dir <- file.path(mcmc.set$meta$output.dir, 'diagnostics')
	if(!file.exists(save.dir)) 
		dir.create(save.dir, recursive=TRUE)
	save.file <- do.call(paste('store.', class.name, sep=''), list(diag, thin=thin, burnin=burnin, 
							output.dir=save.dir))
	if(verbose) cat('\nConvergence diagnostics stored in', save.file, '\n')
	return(diag)
}

tfr.diagnose <- function(sim.dir, thin=80, burnin=2000, express=FALSE, 
							country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
	invisible(.do.diagnose(type='tfr', class.name='bayesTFR.convergence', 
							sim.dir=sim.dir, thin=thin, burnin=burnin, express=express,
							country.sampling.prop=country.sampling.prop, keep.thin.mcmc=keep.thin.mcmc,	verbose=verbose))
}

tfr3.raftery.diag <- function(mcmc=NULL, 
							 sim.dir=file.path(getwd(), 'bayesTFR.output'),
							 burnin=0, country=NULL,
							 par.names = tfr3.parameter.names(),
							 par.names.cs = tfr3.parameter.names.cs(),
							 country.sampling.prop=1,
							 verbose=TRUE, ...) {
	mcmc.set <- if (is.null(mcmc)) get.tfr3.mcmc(sim.dir=sim.dir, low.memory=TRUE) else mcmc
	return(tfr.raftery.diag(mcmc=mcmc.set, burnin=burnin,
						country=country, par.names=par.names, par.names.cs=par.names.cs,
						country.sampling.prop=country.sampling.prop, verbose=verbose, ...))
}

tfr3.diagnose <- function(sim.dir, thin=60, burnin=10000, express=TRUE, 
						country.sampling.prop=NULL, verbose=TRUE, ...) {
	invisible(.do.diagnose(type='tfr3', class.name='bayesTFR.convergence', 
							sim.dir=sim.dir, thin=thin, burnin=burnin, express=express,
							country.sampling.prop=country.sampling.prop, keep.thin.mcmc=FALSE, verbose=verbose))
}

diag.thin.indep <- function(mcmc.list, q) {
	k <- matrix(NA, nrow=length(mcmc.list), ncol=ncol(mcmc.list[[1]]))
	for(imcmc in 1:length(mcmc.list)) {
		mcmc <- mcmc.list[[imcmc]]
		for (i in 1:ncol(mcmc)) 
			k[imcmc,i] <- thinindep(mcmc[,i], q)	}
	return(k)
}

thinindep <- function(x,q){
	## find the smallest integer k that makes the thinned chain independent ## 
	## x is the MCMC samples, and q is the quantile being estimated ##
	k=0
	bic=0
	n0=length(x)
	qx=quantile(x,q)
	while(bic>=0){
		k=k+1
		u=x[seq(1,n0,by=k)]
		n=length(u)
		## change the continuous chain to 0-1 chain ##
		zt <- factor(u<=qx, levels=c(FALSE, TRUE))
		## calculate the 2*2 contingency table ##
		nij <- table(zt[-n], zt[-1])
		nc=apply(nij, 2, sum)
		nr=apply(nij, 1, sum)
		ni=c(nc[1],nc[1],nc[2],nc[2])
		nj=c(nr, nr)
		nij <- c(nij)
		logn <- log(n)
		idx <- nij > 0
		lognij <- log(nij[idx])
		logni <- log(ni[idx])
		lognj <- log(nj[idx])
		## use BIC to determine whether thinning every kth sample is enough ##
		bic <- sum(2*(nij[idx]*(logn+lognij-logni-lognj))) - logn
	}
return(k)
}

tfr.dl.coverage <- function(sim.dir, pi=c(80,90,95), burnin=2000, verbose=TRUE) {
	if(has.tfr.prediction(sim.dir=sim.dir)) {
		pred <- get.tfr.prediction(sim.dir=sim.dir)
		mcmc.set <- pred$mcmc.set
		burnin = 0 # because the prediction mcmc.set is already burned and collapsed
	} else mcmc.set <- get.tfr.mcmc(sim.dir)
	return(.doGoF.dl(mcmc.set, pi=pi, burnin=burnin, verbose=verbose))
}

tfr.DLisDecrement <- function() {
	return(TRUE)
}

.doGoF.dl <- function(mcmc.set, pi=c(80,90,95), type='tfr', burnin=0, verbose=TRUE) {
	meta <- mcmc.set$meta
	countries.index <- get.countries.index(meta)
	data <- get.data.matrix(meta)
	T.total <- nrow(data)
	al.low <- (1 - pi/100)/2
	al.high <- 1 - al.low
	total.GoF <- rep(0, length(pi))
	time.GoF <- matrix(0, nrow=length(pi), ncol=T.total-1)
	country.GoF <- matrix(0, nrow=length(pi), ncol=max(countries.index))
	total.mse <- total.mae <- 0
	time.mse <- time.mae <- rep(0, T.total-1)
	country.mse <- country.mae <- rep(0, max(countries.index))
	counter <- matrix(0, ncol=max(countries.index), nrow=T.total-1)
	pred.cdf <- matrix(NA, ncol=max(countries.index), nrow=T.total-1)
	if(verbose) cat('\nAssessing goodness of fit for country ')
	for(icountry in countries.index) {
		if(verbose) cat(icountry, ', ')
		country.code <- meta$regions$country_code[icountry]
		observed <- diff(data[1:T.total, icountry]) * if(do.call(paste(type,'.DLisDecrement', sep=''), list())) -1 else 1 
		x <- data[1:(T.total - 1), icountry]
		valid.time <- which(!is.na(observed))
		if(length(valid.time) == 0) next
		x <- x[valid.time]
		observed <- observed[valid.time]
		dlc <- do.call(paste(type,'.get.dlcurves', sep=''), 
					list(x, mcmc.set$mcmc.list, country.code, icountry, burnin=burnin, 
							nr.curves=2000, predictive.distr=TRUE))
		counter[valid.time,icountry] <- counter[valid.time,icountry] + 1
		for (i in 1:length(pi)) {
        	dlpi <- apply(dlc, 2, quantile, c(al.low[i], al.high[i]))
        	country.GoF[i,icountry] <- sum(observed >= dlpi[1,] & observed <= dlpi[2,])
        	total.GoF[i] <- total.GoF[i] + country.GoF[i,icountry]
        	for(itime in 1:length(valid.time)) {
        		time.GoF[i,valid.time[itime]] <- time.GoF[i,valid.time[itime]] + (observed[itime] >= dlpi[1,itime] & observed[itime] <= dlpi[2,itime])
        	}
        }
        for(itime in 1:length(valid.time)) {
        	rankdistr <- rank(c(dlc[,itime], observed[itime]))
        	pred.cdf[valid.time[itime], icountry] <- rankdistr[nrow(dlc)+1]/(nrow(dlc)+1)
		}
        dlmean <- apply(dlc, 2, mean)
        dlmedian <- apply(dlc, 2, median)
        country.mse[icountry] <- sum((observed-dlmean)^2)
        country.mae[icountry] <- sum(abs(observed-dlmedian))
        total.mse <- total.mse + country.mse[icountry]
        total.mae <- total.mae + country.mae[icountry]
        for(itime in 1:length(valid.time)) {
        	time.mse[valid.time[itime]] <- time.mse[valid.time[itime]] + (observed[itime] - dlmean[itime])^2
        	time.mae[valid.time[itime]] <- time.mae[valid.time[itime]] + abs(observed[itime] - dlmedian[itime])
        }
	}
	if(verbose) cat('\n')	
	total.GoF <- total.GoF/sum(counter)
	total.mse <- total.mse/sum(counter)
	total.mae <- total.mae/sum(counter)
	pi.names <- paste(pi, '%', sep='')
	names(total.GoF) <- pi.names
	names(total.mse) <- 'RMSE'
	names(total.mae) <- 'MAE'
	rowsum.counter <- rowSums(counter)
	for(row in 1:nrow(time.GoF)) {
		time.GoF[row,] <- time.GoF[row,]/rowsum.counter
	}
	time.mse <- time.mse/rowsum.counter
	time.mae <- time.mae/rowsum.counter
	dimnames(time.GoF) <- list(pi.names, rownames(data)[2:T.total])
	names(time.mse) <- rownames(data)[2:T.total]
	names(time.mae) <- rownames(data)[2:T.total]
	colsum.counter <- colSums(counter)
	for(row in 1:nrow(country.GoF)) {
		country.GoF[row,] <- country.GoF[row,]/colsum.counter
	}
	country.mse <- country.mse/colsum.counter
	country.mae <- country.mae/colsum.counter
	dimnames(country.GoF) <- list(pi.names, meta$regions$country_code[1:max(countries.index)])
	names(country.mse) <- meta$regions$country_code[1:max(countries.index)]
	names(country.mae) <- meta$regions$country_code[1:max(countries.index)]
	country.GoF[is.nan(country.GoF)] <- NA
	country.mse[is.nan(country.mse)] <- NA
	country.mae[is.nan(country.mae)] <- NA
	if(verbose) cat('Done.\n')
	return(list(total.coverage=total.GoF, time.coverage=time.GoF, country.coverage=country.GoF,
				total.rmse=sqrt(total.mse), time.rmse=sqrt(time.mse), country.rmse=sqrt(country.mse), 
				total.mae=total.mae, time.mae=time.mae, country.mae=country.mae,
				pred.cdf=pred.cdf, n=counter))
}