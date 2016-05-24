e0.raftery.diag <- function(mcmc=NULL, 
							 sim.dir=file.path(getwd(), 'bayesLife.output'),
							 burnin=0, country=NULL,
							 par.names = e0.parameter.names(),
							 par.names.cs = e0.parameter.names.cs(),
							 country.sampling.prop=1,
							 verbose=TRUE, ...) {
return(bayesTFR::tfr.raftery.diag(mcmc=mcmc, sim.dir=sim.dir, burnin=burnin,
						country=country, par.names=par.names, par.names.cs=par.names.cs,
						country.sampling.prop=country.sampling.prop, verbose=verbose, ...))
}

e0.diagnose <- function(sim.dir, thin=225, burnin=10000, express=FALSE, 
						country.sampling.prop=NULL, keep.thin.mcmc=FALSE, verbose=TRUE) {
	invisible(bayesTFR:::.do.diagnose(type='e0', class.name='bayesLife.convergence', 
							sim.dir=sim.dir, thin=thin, burnin=burnin, express=express,
							country.sampling.prop=country.sampling.prop, keep.thin.mcmc=keep.thin.mcmc,							verbose=verbose))
	
	
}

e0.dl.coverage <- function(sim.dir, pi=c(80,90,95), burnin=10000, verbose=TRUE) {
	if(has.e0.prediction(sim.dir=sim.dir)) {
		pred <- get.e0.prediction(sim.dir=sim.dir)
		mcmc.set <- pred$mcmc.set
		burnin = 0 # because the prediction mcmc.set is already burned and collapsed
	} else mcmc.set <- get.e0.mcmc(sim.dir)
	return(.do.e0GoF.dl(mcmc.set, pi=pi, burnin=burnin, verbose=verbose))
}

e0.DLisDecrement <- function() {
	return(FALSE)
}

.do.e0GoF.dl <- function(mcmc.set, pi=c(80,90,95), burnin=0, verbose=TRUE) {
	# Like .doGoF.dl (in bayesTFR) but it also considers supplemental data
	meta <- mcmc.set$meta
	countries.index <- get.countries.index(meta)
	data <- get.data.matrix(meta)
	T.total <- nrow(data)
	T.recent.start <- 0
	T.names <- rownames(data)[1:(T.total-1)]
	if(!is.null(meta$suppl.data$e0.matrix)) {
		T.recent.start <- nrow(meta$suppl.data$e0.matrix) 
		T.total <- T.total + T.recent.start
		T.names <- c(rownames(meta$suppl.data$e0.matrix), T.names)
	}
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
		valid.time <- which(!is.na(meta$d.ct[,icountry]))
		if(length(valid.time) == 0) next
		observed <- meta$d.ct[valid.time, icountry]
		x <- data[valid.time, icountry]
		valid.time.all <- valid.time + T.recent.start
		if(!is.null(meta$suppl.data$e0.matrix)) {
    		supp.c.idx <- which(is.element(meta$suppl.data$index.to.all.countries, icountry))
    		if(length(supp.c.idx) > 0) {
    			suppl.data.idx <- which(!is.na(meta$suppl.data$d.ct[,supp.c.idx]))
    			x <- c(meta$suppl.data$e0.matrix[suppl.data.idx, supp.c.idx], x)
    			observed <- c(meta$suppl.data$d.ct[suppl.data.idx, supp.c.idx], observed)
    			valid.time.all <- c(suppl.data.idx, valid.time.all)
   			}
    	}
		dlc <- e0.get.dlcurves(x, mcmc.set$mcmc.list, country.code, icountry, burnin=burnin, 
							nr.curves=2000, predictive.distr=TRUE)
		counter[valid.time.all,icountry] <- counter[valid.time.all,icountry] + 1
		for (i in 1:length(pi)) {
        	dlpi <- apply(dlc, 2, quantile, c(al.low[i], al.high[i]))
        	country.GoF[i,icountry] <- sum(observed >= dlpi[1,] & observed <= dlpi[2,])
        	total.GoF[i] <- total.GoF[i] + country.GoF[i,icountry]
        	for(itime in 1:length(valid.time.all)) {
        		time.GoF[i,valid.time.all[itime]] <- time.GoF[i,valid.time.all[itime]] + (observed[itime] >= dlpi[1,itime] & observed[itime] <= dlpi[2,itime])
        	}
        }
        for(itime in 1:length(valid.time.all)) {
        	rankdistr <- rank(c(dlc[,itime], observed[itime]))
        	pred.cdf[valid.time.all[itime], icountry] <- rankdistr[nrow(dlc)+1]/(nrow(dlc)+1)
		}
        dlmean <- apply(dlc, 2, mean)
        dlmedian <- apply(dlc, 2, median)
        country.mse[icountry] <- sum((observed-dlmean)^2)
        country.mae[icountry] <- sum(abs(observed-dlmedian))
        total.mse <- total.mse + country.mse[icountry]
        total.mae <- total.mae + country.mae[icountry]
        for(itime in 1:length(valid.time.all)) {
        	time.mse[valid.time.all[itime]] <- time.mse[valid.time.all[itime]] + (observed[itime] - dlmean[itime])^2
        	time.mae[valid.time.all[itime]] <- time.mae[valid.time.all[itime]] + abs(observed[itime] - dlmedian[itime])
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
	dimnames(time.GoF) <- list(pi.names, T.names)
	names(time.mse) <- T.names
	names(time.mae) <- T.names
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
