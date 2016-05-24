

runjags.dic <- function(deviance.table, deviance.sum, mcmclist){
	
	# If deviance in mcmclist use it for chains
	
	if(runjags.getOption('debug')>=10){
		cat('Deviance table and deviance sum vector:\n')
		print(deviance.table)
		print(deviance.sum)
	}
	
	if(all(is.na(deviance.table)) && all(is.na(deviance.sum))){
		dic <- "DIC statistics not available"		
	}else{
		
		if('deviance' %in% varnames(mcmclist)){
			meandeviance.chain <- sapply(mcmclist, function(x) return(mean(x[,'deviance'])))			
		}else{
			meandeviance.chain <- rep(NA, nchain(mcmclist))	
		}

    	meandeviance <- deviance.sum['sum.mean.deviance']
		meanpopt <- deviance.sum['sum.mean.pOpt']
		meanpd <- deviance.sum['sum.mean.pD']
		
		dicstats <- list(dic=meandeviance+meanpd, dic.chains=meandeviance.chain+meanpd, ped=meandeviance+meanpopt, ped.chains=meandeviance.chain+meanpopt, meandeviance=meandeviance, meandeviance.chains=meandeviance.chain, meanpd=meanpd, meanpopt=meanpopt)

		class(dicstats) <- 'dicstats'
		dic <- dicstats
	}
	
}



runjags.summaries <- function(fullmcmclist, thinnedmcmclist, psrf.target, normalise.mcmc, modeest.opts, confidence, autocorr.lags, custom, silent=FALSE){
	
	n.chains <- coda::nchain(thinnedmcmclist)
	n.iter <- coda::niter(thinnedmcmclist)
	n.var <- coda::nvar(thinnedmcmclist)
	
	# 1 iteration + 1 chain is handled by normalise.mcmcfun:
	normalised <- normalise.mcmcfun(thinnedmcmclist, normalise = normalise.mcmc, warn=TRUE, remove.nonstochastic = TRUE)
	
	normalisedmcmc <- normalised$mcmc
	truestochastic <- normalised$truestochastic
	semistochastic <- normalised$semistochastic
	nonstochastic <- normalised$nonstochastic			
	
	collapsed <- combine.mcmc(thinnedmcmclist, collapse.chains=TRUE)
	
	options(show.error.messages = FALSE)
	success <- try({
	  suppressWarnings(tsummary <- summary(combine.mcmc(thinnedmcmclist, collapse.chains=FALSE)))
	  if(class(tsummary$statistics)=="numeric"){
	    tsummary$statistics <- t(as.matrix(tsummary$statistics))
	    dimnames(tsummary$statistics)[[1]] <- varnames(thinnedmcmclist)
	    tsummary$quantiles <- t(as.matrix(tsummary$quantiles))
	    dimnames(tsummary$quantiles)[[1]] <- varnames(thinnedmcmclist)
	  }
	})
	if(inherits(success, 'try-error')) tsummary <- "An unexpected error occured while calculating summary statistics"
	options(show.error.messages = TRUE)			
	
	# First lot of summaries require at least one stochastic variable:
	if(!all(nonstochastic)){
	
		autocorr.lags <- autocorr.lags[autocorr.lags < n.iter]
		if(length(autocorr.lags)==0){
			autocorr.lags <- 1
		}
	
		autocorrelation <- safe.autocorr.diag(fullmcmclist[,!nonstochastic,drop=FALSE], lags=autocorr.lags)
		suppressWarnings(crosscorrelation <- crosscorr(fullmcmclist[,!nonstochastic,drop=FALSE]))
		class(crosscorrelation) <- "crosscorrstats"	
	
		success <- try({
		
			if(n.chains > 1 && n.iter > 1){
				
				if(!silent) swcat("Calculating the Gelman-Rubin statistic for ", nvar(thinnedmcmclist), " variables....\n", sep="")
				convergence <- safe.gelman.diag(normalisedmcmc, transform=FALSE, autoburnin=FALSE)
		
				convergence <- c(convergence, psrf.target=psrf.target)
				class(convergence) <- "gelmanwithtarget"
		
				#n.params <- nrow(convergence$psrf)
				
				if(nrow(convergence$psrf) != sum(!nonstochastic))
					stop(paste(nrow(convergence$psrf), ' statistics were returned by gelman.diag but ', sum(!nonstochastic), ' were expected', sep=''))

			}else{
				if(n.iter==1){
					if(runjags.getOption('summary.warning'))
						warning("Convergence cannot be assessed with only 1 iteration", call.=FALSE)
					
					convergence <- "Convergence cannot be assessed using only 1 iteration"
				}
				if(n.chains==1){
					if(runjags.getOption('summary.warning'))
						warning("Convergence cannot be assessed with only 1 chain", call.=FALSE)
					
					convergence <- "Convergence cannot be assessed using only 1 chain"
				}
				
				#param.conv <- 1
				#n.params <- 1
			}

			if(n.chains > 1 && n.iter > 1){
				if(class(convergence$mpsrf)!="numeric"){
					mpsrfstring <- " (Unable to calculate the multi-variate psrf)"
				}else{
					mpsrfstring <- paste(" (multi-variate psrf = ", round(convergence$mpsrf, digits=3), ")", sep="")
				}
			}
		
			##########################################################
			#### REMOVED CODE
			##########################################################		
	#		autocorrelated <- 0
	#		unconverged <- 0
	#		crosscorrelated <- 0
	#		
	#		for(j in 1:n.params){
	#			if(n.chains > 1){
	#				param.conv <- convergence$psrf[j, 1]
	#				if(!is.na(param.conv)){
	#					if(param.conv > psrf.target){
	#						unconverged <- unconverged + 1
	#					}
	#				}
	#			}
	#			param.autocorr <- autocorrelation[3,j]
	#			if(!is.na(param.autocorr)){
	#				if(param.autocorr > 0.1){
	#					autocorrelated <- autocorrelated + 1
	#				}
	#			}
	#			param.crosscorr <- crosscorrelation
	#			param.crosscorr[1:nrow(param.crosscorr), 1:ncol(param.crosscorr)] <- 0
	#			# print("0.3 is an arbitrary figure for crosscorr")
	#			# param.crosscorr is symmetrical - so divide by 2 to get number cross correlated:
	#			crosscorrelated <- sum(param.crosscorr > 0.3)/2
	#		}
	#
	#		updates <- niter(mcmclist)
	#
	#		if(!is.na(param.conv)){
	#			if(unconverged > 0){
	#				if(n.params==1 & !silent) swcat("Convergence may have failed for this run after ", updates, " iterations (psrf = ", round(convergence$psrf[1,1], digits=3), ")\n", sep="") else swcat("Convergence may have failed for this run for ", unconverged, " parameter", if(unconverged>1) "s", " after ", updates, " iterations", mpsrfstring, "\n", sep="")
	#			}else{
	#				if(n.chains > 1 & !silent) swcat("The Gelman-Rubin statistic is below ", psrf.target, " for all parameters\n", sep="")
	#				if(n.chains==1 & !silent) swcat("Calculating the Gelman-Rubin statistic requires two or more chains\n")
	#			}
	#		}else{
	#			if(!silent) swcat("There was an unexpected error calculating the Gelman-Rubin statistic\n")
	#		}
	#
	#		if(!is.na(param.autocorr)){
	#			if(autocorrelated > 0){
	##				if(!silent) swcat("IMPORTANT:  There was a high degree of autocorrelation for ", autocorrelated, " parameter", if(autocorrelated>1) "s", " (see the $autocorr element of the runjags object for more details)\n", sep="")
	#			}
	#		}else{
	#			if(!silent) swcat("There was an unexpected error calculating the autocorrelation dependence\n")
	#		}
	#		if(crosscorrelated > 0){
	##			if(!silent) swcat("IMPORTANT:  There was a high degree of cross-correlation for ", crosscorrelated, " parameter pair", if(crosscorrelated>1) "s", " (see the $crosscorr element of the runjags object for more details)\n", sep="")
	#		}
		##########################################################
		##########################################################
	
	
		}, silent=FALSE)
	
		if(inherits(success, 'try-error')){
			if(runjags.getOption('debug'))
				stop("An unexpected error occured when assessing convergence")
			if(!silent)
				swcat("An unexpected error occured when assessing convergence\n")
			convergence <- "An unexpected error occured when assessing convergence"
		}
	

		s <- try(sseff <- effectiveSize(thinnedmcmclist), silent=TRUE)
		if(class(s)=='try-error'){
			if(runjags.getOption('summary.warning'))
				warning('There was an error calculating the effective sample size [using coda::effectiveSize()] for one or more parameters', call.=FALSE)
			
			sseff <- apply(collapsed,2,function(x){
				ess <- try(size <- effectiveSize(x))
				if(class(ess)=='try-error')
					return(NA)
				else
					return(size)
			})
		}
	
		if(any(confidence>1 | confidence<0)){
			runjags.getOption('summary.warning')
				warning('Invalid value for confidence was ignored (this must be in the range 0-1)')
			confidence <- confidence[confidence<1 & confidence>0]
			if(length(confidence)==0) confidence <- 0.95
		}
		confidence <- sort(confidence, decreasing=TRUE)
		nc <- length(confidence)
	
		options(show.error.messages = FALSE)
		success <- try({
			thpd <- matrix(0, ncol=1+(length(confidence)*2), nrow=n.var)
			thpd[,nc+1] <- apply(collapsed,2,median)
			for(i in 1:nc){
				suppressWarnings(thpd[,c(i, (nc*2 +2)-i)] <- HPDinterval(collapsed, prob=confidence[i]))
			}
			dimnames(thpd) <- list(varnames(collapsed), c(paste("Lower",round(confidence*100),sep=""), "Median", paste("Upper",round(confidence[nc:1]*100),sep="")))
		})

		if(inherits(success, 'try-error')) thpd <- "An unexpected error occured while calculating summary statistics"
		options(show.error.messages = TRUE)			

		options(show.error.messages = FALSE)

		success <- try({
	#		stochastic <- tsummary$statistics[,2] != 0
	#		sseff <- effectiveSize(collapsed)
			se <- tsummary$statistics[,2]
			mcse <- se / sqrt(sseff)
			sseff <- sseff[!nonstochastic]
			mcse <- mcse[!nonstochastic]
			se <- se[!nonstochastic]
			thmcse <- list(sseff=sseff, ssd=se, mcse=mcse)
			})	
		if(inherits(success, 'try-error')) thmcse <- "An unexpected error occured while calculating Monte Carlo error"
		options(show.error.messages = TRUE)			
	
		class(thmcse) <- 'mcsestats'
		
		
	}else{
		
		autocorrelation <- matrix(NA, ncol=n.var, nrow=1, dimnames=list('Lag.10', NULL))
		crosscorrelation <- NA
		convergence <- list()
		sseff <- NA
		confidence <- NA
		thpd <- matrix(rep(tsummary$statistics[,1], 2), ncol=2, nrow=n.var, dimnames=list(NULL, c('Lower95','Upper95')))
		thmcse <- list(sseff=NA, ssd=NA, mcse=NA)
		
	}


	# Calculate the mode:
	modestats <- numeric(n.var)
	success <- try({
	modestats <- apply(collapsed,2,function(x){
	  if(all(abs(x - round(x, digits=0)) < .Machine$double.eps^0.5)){
	    return(as.numeric(names(sort(table(x), decreasing=TRUE))[1]))
	  }else{
	    return(NA)
	  }
	})
	
	discrete <- !is.na(modestats)
	
	if(any(is.na(modestats)) && runjags.getOption('mode.continuous')){
		if(!suppressPackageStartupMessages(requireNamespace('modeest', quietly=TRUE)))
			stop('The "modeest" package is required to calculate the mode of continuous variables', call.=FALSE)
		if(is.null(modeest.opts)) modeest.opts <- list()
		if(!is.list(modeest.opts)){
			warning('Non list value provided for modeest.opts was ignored')
			modeest.opts <- list()
		}
		if(!any(names(modeest.opts)=='method'))
			modeest.opts$method <- 'shorth'
		if(!any(names(modeest.opts)=='tie.limit'))
			modeest.opts$tie.limit <- Inf
		modeest.opts$MARGIN <- 2
		modeest.opts$FUN <- function(x, ...){
			return(modeest::mlv(x, ...)$M[1])
		}
		modeest.opts$X <- collapsed[,is.na(modestats)&!nonstochastic,drop=FALSE]
		modestats[is.na(modestats)&!nonstochastic] <- do.call('apply', args=modeest.opts)
	}
	})
	if(inherits(success, 'try-error')) warning('An unexpected error occured while calculating the mode')
	
	
	# Possible custom function:
	customstats <- matrix(NA, ncol=1, nrow=n.var)
	if(is.function(custom)){
		success <- try({
			customstats <- apply(collapsed,2,custom)
			if(is.null(dim(customstats))){
				customstats <- matrix(customstats,nrow=1,dimnames=list('Custom.Summary',names(customstats)))
			}
			if(length(dim(customstats)!=2) && dim(customstats)[2]!=n.var){
				warning('Custom summary produced an incompatible result - each variable must have the same number of summary statistics returned')
				customstats <- matrix(NA, ncol=1, nrow=n.var)
			}
			if(is.null(dimnames(customstats)[[1]])) dimnames(customstats) <- list(paste('Custom.',1:nrow(customstats),sep=''), dimnames(customstats)[[2]])
			customstats <- t(customstats)
		})
		
		if(inherits(success, 'try-error')) warning('An unexpected error occured while calculating the custom summary function')
	}	
	
	# Pre-create summary table:
	tsummaries <- cbind(thpd, tsummary$statistics[,1:2,drop=FALSE], modestats)
  	sumnames <- c(dimnames(thpd)[[2]], "Mean","SD","Mode")

	if(is.function(custom)){
	tsummaries <- cbind(tsummaries, customstats)
	sumnames <- c(sumnames, dimnames(customstats)[[2]])
	}  

	# Add 'NA' for non-stochastic variables:
	mcse=pos=sse=psrfs=psrfs.upper <- replicate(length(nonstochastic),NA)
	mcse[!nonstochastic] <- thmcse$mcse
	sse[!nonstochastic] <- round(thmcse$sse)
	pos[!nonstochastic] <- round((thmcse$mcse/tsummary$statistics[!nonstochastic,2])*100,1)
	
	autocorrs <- matrix(ncol=nrow(autocorrelation), nrow=length(!nonstochastic))
	autocorrs[!nonstochastic,] <- t(autocorrelation)
	
	# Catch if only 1 chain | iteration
	if(n.chains==1 || n.iter==1){
	  psrfs <- replicate(length(nonstochastic), NA)
	  psrfs.upper <- replicate(length(nonstochastic), NA)  
	} else {
	  psrfs[!nonstochastic] <- convergence$psrf[,1]
	  psrfs.upper[!nonstochastic] <- convergence$psrf[,2]
	}

	tsummaries <- cbind(tsummaries, mcse, pos, sse, autocorrs, psrfs)#, psrfs.upper)
	sumnames <- c(sumnames, "MCerr","MC%ofSD","SSeff",gsub("Lag ", "AC.", dimnames(autocorrelation)[[1]]), "psrf")#, "psrf.UCI")
	
	dimnames(tsummaries) <- list(dimnames(tsummaries)[[1]], sumnames)
	
	return(list(summaries=tsummaries, summary=tsummary, HPD=thpd, hpd=thpd, mcse=thmcse, psrf=convergence, autocorr=autocorrelation, crosscorr=crosscorrelation, truestochastic=truestochastic, semistochastic=semistochastic, nonstochastic=!(truestochastic | semistochastic), discrete=discrete))	
	
}


runjagsplots <- function(mcmclist, psrfs, discrete, silent=FALSE, trace=TRUE, density=TRUE, histogram=TRUE, ecdf=TRUE, autocorr=TRUE, crosscorr=TRUE, key=TRUE, col=NA, trace.iters=1000, separate.chains=FALSE, trace.options=list(), density.options=list(), histogram.options=list(), ecdfplot.options=list(), acplot.options=list()){
	
	if(niter(mcmclist) < 2)
		stop('Unable to produce plots from less than 2 iterations', call.=FALSE)
	
	# Implement separate.chains at some point.... as an element of trace/density/ecdf option lists
	# separate.chains <- FALSE
	
	# Use default colours from lattice:
	cols <- c("#0080ff", "#ff00ff", "darkgreen", "#ff0000", "orange", "#00ff00", "brown")
	if(nchain(mcmclist)>length(cols)) cols <- rainbow(nchain(mcmclist)+1)
	cols <- cols[1:nchain(mcmclist)]
	cols <- c(cols, 'dark grey')
	
	if(!identical(col, NA)){
		if(length(col)==nchain(mcmclist)){
			col <- c(col, 'dark grey')
		}
		cols <- numeric(nchain(mcmclist)+1)
		cols[] <- col
	}
	
	if(!length(discrete)==nvar(mcmclist))
		stop('An error occured while creating plots: the length of the discrete identifier did not match the number of variables - please file a bug report to the package author')
	
	if(any(names(trace.options)=='col')) warning('The "col" argument in trace.options was ignored')
	
	success <- try({

	if(trace || density || ecdf || histogram || autocorr){
		
		thinned.mcmc <- combine.mcmc(list(mcmclist), collapse.chains=FALSE, return.samples=min(trace.iters, niter(mcmclist)))
		if(separate.chains){
			nplots <- nchain(thinned.mcmc)
			chainnames <- paste('chain_', gsub(' ', '0', format(1:nplots, scientific=FALSE)), sep='')
		}else{
			nplots <- 1
			chainnames <- 'chains_all'
		} 
		niters <- niter(thinned.mcmc)
		
		plot1 = plot2 = plot3 = plot4 = plot5 <- lapply(1:length(varnames(thinned.mcmc)), function(x){
			newl <- vector('list', length=nplots)
			names(newl) <- chainnames
			return(newl)
			}) 
		names(plot1) = names(plot2) = names(plot3) = names(plot4) = names(plot5)  <- varnames(thinned.mcmc)
		
		# Was using this to manually over-write the x axis labels to start at burnin+1, but it's too complicated:
#		iternames <- dimnames(thinned.mcmc[[1]])[[1]]
#		iterrange <- as.numeric(c(iternames[1],iternames[length(iternames)]))
#		ilabels <- signif(seq(iterrange[1], iterrange[2], length.out=5)+1,3)
#		iat <- seq(0,niter(thinned.mcmc),length.out=5)
		
		for(i in 1:length(varnames(thinned.mcmc))){
			# Too flowery:  paste("Value of '", dimnames(plotdata[[1]])[[2]], "'", sep="")
			all.thinned <- as.mcmc.list(lapply(thinned.mcmc, function(x) return(x[,i,drop=FALSE]))) # xyplot throws an error if not a matrix
			all.full <- as.mcmc.list(lapply(mcmclist, function(x) return(x[,i,drop=FALSE]))) # xyplot throws an error if not a matrix			
			varname <- dimnames(all.thinned[[1]])[[2]]
			
			for(p in 1:nplots){				
				if(separate.chains){
					thinplotdata <- all.thinned[[p]]
					allplotdata <- all.full[[p]]
					chain <- p
					nchains <- 1
					varlab <- paste(varname, ' (chain ', chain, ')', sep='')
					usecols <- cols[p]
				}else{
					thinplotdata <- all.thinned
					allplotdata <- all.full
					chain <- paste(p,':',nplots,sep="")
					varlab <- varname
					nchains <- nchain(thinplotdata)
					usecols <- cols[1:nchains]
				}
				combcol <- cols[nchains+1]
				
				if(trace){
				
					plotopts <- trace.options
					plotopts$col <- usecols
					plotopts$varname <- varname
				
					if(!any(names(plotopts)=='ylab')) plotopts$ylab <- varlab
					if(!any(names(plotopts)=='xlab')) plotopts$xlab <- "Iteration"
				
					plotopts$ylab <- eval(plotopts$ylab)
					plotopts$xlab <- eval(plotopts$xlab)
				
					if(any(names(plotopts)=='main')) plotopts$main <- eval(plotopts$main)
					if(any(names(plotopts)=='sub')) plotopts$sub <- eval(plotopts$sub)
				
					plotopts$x <- thinplotdata
					plot1[[i]][[p]] <- do.call('xyplot', args=plotopts)
					# plot1[[i]] <- xyplot(plotdata, ylab=dimnames(plotdata[[1]])[[2]], xlab="Iteration", col=cols)  #doesn't work: , scales=list(x=list(at=iat, labels=ilabels)))
				
				}
				if(density){
					plotopts <- density.options
					plotopts$col <- usecols
					plotopts$varname <- varname
				
					if(!any(names(plotopts)=='ylab')) plotopts$ylab <- "Density"
					if(!any(names(plotopts)=='xlab')) plotopts$xlab <- varlab
					if(!any(names(plotopts)=='aspect')) plotopts$aspect <- 'fill'
					if(!any(names(plotopts)=='plot.points')) plotopts$plot.points <- FALSE
				
					plotopts$ylab <- eval(plotopts$ylab)
					plotopts$xlab <- eval(plotopts$xlab)
				
					if(any(names(plotopts)=='main')) plotopts$main <- eval(plotopts$main)
					if(any(names(plotopts)=='sub')) plotopts$sub <- eval(plotopts$sub)
				
					plotopts$x <- allplotdata
					plot2[[i]][[p]] <- do.call('densityplot', args=plotopts)
					# plot2[[i]] <- densityplot(plotdata, plot.points=FALSE, ylab="Density", xlab=dimnames(plotdata[[1]])[[2]], col=cols, aspect='fill')
				} 
				if(ecdf){		
					plotopts <- ecdfplot.options
					plotopts$col <- usecols
					plotopts$varname <- varname
				
					if(!any(names(plotopts)=='ylab')) plotopts$ylab <- "ECDF"
					if(!any(names(plotopts)=='xlab')) plotopts$xlab <- varlab
				
					plotopts$ylab <- eval(plotopts$ylab)
					plotopts$xlab <- eval(plotopts$xlab)
				
					if(any(names(plotopts)=='main')) plotopts$main <- eval(plotopts$main)
					if(any(names(plotopts)=='sub')) plotopts$sub <- eval(plotopts$sub)

						
					# niters is what to thin to:
					tdat <- data.frame(cd=1:niters/niters, vx=numeric(niters * max(1, nchains*!separate.chains)))
					if(separate.chains){
						tdat$vx <- quantile(as.numeric(unlist(allplotdata)), probs=tdat$cd)
					}else{
						chains <- rep(1:nchains, each=niters)
						for(c in 1:nchains) tdat$vx[chains==c] <- as.numeric(quantile(as.numeric(allplotdata[[c]]), probs=1:niters/niters))
						plotopts$groups <- chains
					}
						
					plotopts$data <- tdat
					plotopts$x <- cd ~ vx
					
					if(!any(names(plotopts)=='type')) plotopts$type <- 's'
					if(!any(names(plotopts)=='ylim')) plotopts$ylim <- c(-0.05,1.05)
					if(!any(names(plotopts)=='panel')) plotopts$panel <- function(x,y,...){
						panel.abline(h=0)	
						panel.abline(h=1)
						panel.xyplot(x,y,...)
					}
					
					plot3[[i]][[p]] <- do.call('xyplot', args=plotopts)
				} 
				
				# Either way the rest are all the full data (not thinned):
				if(!separate.chains){
					stopifnot(p==1)
					plotdata <- combine.mcmc(allplotdata, collapse.chains=TRUE)
				}else{
					plotdata <- allplotdata
				}

				if(histogram){		
					plotopts <- histogram.options
					plotopts$varname <- varname

					if(separate.chains){
						plotopts$border <- usecols
						plotopts$col <- usecols
					}
					if(!any(names(plotopts)=='border')) plotopts$border <- combcol
					if(!any(names(plotopts)=='col')) plotopts$col <- combcol
				
					if(!any(names(plotopts)=='xlab')) plotopts$xlab <- varlab
					if(!any(names(plotopts)=='ylab')) plotopts$ylab <- "% of total"
				
					plotopts$ylab <- eval(plotopts$ylab)
					plotopts$xlab <- eval(plotopts$xlab)
				
					if(any(names(plotopts)=='main')) plotopts$main <- eval(plotopts$main)
					if(any(names(plotopts)=='sub')) plotopts$sub <- eval(plotopts$sub)
					
					if(!discrete[i] && !any(names(plotopts)=='breaks'))
						plotopts$breaks <- 100
					
					# Maybe use a barchart if < 20 groups or something?
					# Doesn't seem necessary - histogram looks OK for dichotomous at least
#					if(discrete[i] && !is.null(plotopts$force.histogram) && plotopts$force.histogram <= length(unique(plotdata))){
#						tabs <- table(factor(plotdata))
#						plotopts$data <- data.frame(x=as.numeric(names(tabs)), y=tabs)
#						plotopts$x <- y ~ x
#						fcall <- 'barchart'
#					}else{
						plotopts$data <- data.frame(x=as.numeric(unlist(plotdata)))
						plotopts$x <-  ~ x
						fcall <- 'histogram'
#					}
					
					plot4[[i]][[p]] <- do.call(fcall, args=plotopts)
				}
				
				if(autocorr){		
					plotopts <- acplot.options
					plotopts$varname <- varname

					if(separate.chains){
						plotopts$border <- usecols
						plotopts$col <- usecols
					}
					if(!any(names(plotopts)=='border')) plotopts$border <- combcol
					if(!any(names(plotopts)=='col')) plotopts$col <- combcol
				
					if(!any(names(plotopts)=='ylab')) plotopts$ylab <- paste('Autocorrelation of ', varlab, sep='')
					if(!any(names(plotopts)=='xlab')) plotopts$xlab <- 'Lag'
				
					plotopts$ylab <- eval(plotopts$ylab)
					plotopts$xlab <- eval(plotopts$xlab)
				
					if(any(names(plotopts)=='main')) plotopts$main <- eval(plotopts$main)
					if(any(names(plotopts)=='sub')) plotopts$sub <- eval(plotopts$sub)
					
					plotopts$horizontal <- FALSE
					if(!any(names(plotopts)=='origin')) plotopts$origin <- 0
					if(!any(names(plotopts)=='ylim.ratio')) plotopts$ylim <- c(-1.05,1.05)
					if(!any(names(plotopts)=='box.ratio')) plotopts$box.ratio <- 0.5
					
					acfs <- acf(plotdata, plot=FALSE)
					labels <- (0:floor(length(acfs$lag)/5))*5
					if(!any(names(plotopts)=='box.scales')) plotopts$scales <- list(x=list(labels=labels, at=labels+1, tck=1))
					
					plotopts$data <- data.frame(acfs=as.numeric(acfs$acf), lags=as.numeric(acfs$lag))
					plotopts$x <-  acfs ~ lags
					
					plot5[[i]][[p]] <- do.call('barchart', args=plotopts)

				}
			}
		}
		
		class(plot1) <- 'runjagsplots'
		class(plot2) <- 'runjagsplots'
		class(plot3) <- 'runjagsplots'
		class(plot4) <- 'runjagsplots'
		class(plot5) <- 'runjagsplots'
		
		#if(!is.null(startdev)){
		#	for(i in dev.list()){
		#		if(!any(startdev==i)) dev.off(i)
		#	}
		#}else{
		#	try(a <- dev.off(), silent=TRUE)
		#}
		
	}
	
	if(key){
		nchains <- nchain(mcmclist)
		ys <- -(1:(nchains+(!separate.chains && (autocorr || histogram))) -0.5)
		yls <- paste('Chain ', 1:nchains, sep='')
		if((!separate.chains && (autocorr || histogram))) yls <- c(yls, 'Combined')

		kp <- xyplot(0~0, xlim=c(0,1), ylab=NULL, xlab=NULL, scales=list(x=list(labels=NULL, at=NULL), y=list(labels=yls, at=ys)), ylim=c(-(nchains+(!separate.chains && (autocorr || histogram))),0), panel=function(){
		  for(i in 1:nchains){
		    panel.lines(y=-(i-0.5), x=c(0.1,0.9),col=cols[i],lwd=5)
		  }
		  if((!separate.chains && (autocorr || histogram))){
			i <- nchains+1
  		   	panel.lines(y=-(i-0.5), x=c(0.1,0.9),col=cols[i],lwd=5)
		  }
		})
		keyplot <- list(key=kp)
		class(keyplot) <- 'runjagsplots'
	}
	
	if(crosscorr){
		ccplot <- list(crosscorr=crosscorr.plot.lattice(mcmclist))
		class(ccplot) <- 'runjagsplots'
	}		

	})
	if(inherits(success, 'try-error')){
		trace=density=autocorr=crosscorr=key=histogram=ecdf <- FALSE
		warning("An unexpected error occured while attempting to plot graphs")
	}
	
	if(!trace) plot1 <- 'No pre-drawn trace plots available'
	if(!density) plot2 <- 'No pre-drawn density plots available'
	if(!key) keyplot <- 'No pre-drawn density plots available'
	if(!autocorr) plot5 <- 'No pre-drawn autocorrelation plots available'
	if(!crosscorr) ccplot <- 'No pre-drawn crosscorrelation plots available'
	if(!histogram) plot4 <- 'No pre-drawn autocorrelation plots available'
	if(!ecdf) plot3 <- 'No pre-drawn crosscorrelation plots available'

	return(list(trace=plot1, density=plot2, ecdfplot=plot3, histogram=plot4, acplot=plot5, key=keyplot, ccplot=ccplot))
}




# crosscorr.plot function modified from coda version 0.16-1, to switch to lattice based graphics
# Original functions are GPL>=2, copyright of the authors:  Martyn Plummer [aut, cre, trl], Nicky Best [aut], Kate Cowles [aut], Karen Vines [aut], Deepayan Sarkar [aut], Russell Almond [ctb]

crosscorr.plot.lattice <- function (x, col = topo.colors(10)) 
{
    Nvar <- nvar(x)
    pcorr <- crosscorr(x)
    dens <- ((pcorr + 1) * length(col))%/%2 + (pcorr < 1) + (pcorr < -1)
    cutoffs <- format(seq(from = 1, to = -1, length = length(col) + 1), digits = 2)
    leg <- paste("(", cutoffs[-1], ",", cutoffs[-length(cutoffs)], "]", sep = "")
	
    yval <- seq(from = Nvar/2, to = Nvar, length = length(col) + 1)
    ydelta <- Nvar/(2 * (length(col) + 1))

	ay <- list(labels=abbreviate(varnames(x, allow.null=FALSE), minlength=7), at=(1:Nvar)-0.5, rot=90)
	ax <- list(labels=abbreviate(varnames(x, allow.null=FALSE), minlength=7), at=(Nvar:1)-0.5)
	lims <- c(-0.1, Nvar+0.1)
	replot <- xyplot(0 ~ 0, xlim=lims, ylim=lims, xlab = "", ylab = "", scales=list(x=ay, y=ax), panel=function(){
		
		for(i in 1:Nvar){
			ypoints <- c(i-1,i-1,i,i)
			for(j in 1:(Nvar-i+1)){
				panel.polygon(y=ypoints, x=c(j-1,j,j,j-1), col=col[dens[nrow(dens) - i + 1, j]])
			}
		}
		for(i in 1:length(col)){
			panel.polygon(y = c(yval[i], yval[i + 1], yval[i + 1], yval[i], 
            yval[i]), col = col[i], x = c(Nvar - ydelta, Nvar - 
            ydelta, Nvar, Nvar, Nvar - ydelta))
		}
	    panel.text(Nvar - (ydelta+0.07), Nvar, "1", adj = c(1, 1))
	    panel.text(Nvar - (ydelta+0.07), 0.5 * Nvar, "-1", adj = c(1, 0))
	    panel.text(Nvar - (ydelta+0.07), 0.75 * Nvar, "0", adj = c(1, 0.5))
		
	})
	
	return(replot)
}
