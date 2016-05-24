if(getRversion() >= "2.15.1") utils::globalVariables("UNlocations")

pop.aggregate <- function(pop.pred, regions, input.type=c('country', 'region'),
						name = input.type,
						inputs=list(e0F.sim.dir=NULL, e0M.sim.dir='joint_', tfr.sim.dir=NULL),
						my.location.file=NULL,
						verbose=FALSE, ...) {
	if(is.null(my.location.file))
		bayesTFR:::load.bdem.dataset('UNlocations', pop.pred$wpp.year, envir=globalenv(), verbose=verbose)
	else {
		env <- globalenv()
		assign("UNlocations", read.delim(my.location.file, comment.char='#'), envir=env)
	}
	regions <- unique(regions)
	method <- match.arg(input.type)
	if(missing(name)) name <- method
	if(method == 'country') 
		aggr.pred <- pop.aggregate.countries(pop.pred, regions, name, verbose=verbose, ...)
	if(method == 'region')
		aggr.pred <- pop.aggregate.regional(pop.pred, regions, name, inputs=inputs, verbose=verbose)
	invisible(get.pop.aggregation(pop.pred=pop.pred, name=name))
}

get.countries.for.region <- function(region, pop.pred) {
	reg.idx <- which(UNlocations[,'country_code'] == region)
	if(length(reg.idx) <= 0) {
		warning('No region code ', region, ' available.')
		return(c())
	}
	all.countries <- UNlocations[UNlocations[,'location_type'] == 4,]
	location.type <- UNlocations[reg.idx,'location_type']
	if(location.type == 0)  # corresponds to world, i.e. all countries
		countries <- all.countries[,'country_code']
	else {
		code.column <- switch(as.character(location.type), '2'='area_code', '3'='reg_code', paste0('agcode_', location.type))
		if(!is.element(code.column, colnames(UNlocations))) {
			warning('Invalid location type ', location.type, ' for region ', region,'. Location file must contain column ', code.column, '.')
			return(c())
		}	
		countries <- all.countries[is.element(all.countries[,code.column], region), 'country_code']
	}
	return(countries[is.element(countries, pop.pred$countries[,'code'])])
}

pop.aggregate.regional <- function(pop.pred, regions, name,
						inputs=list(e0F.sim.dir=NULL, e0M.sim.dir='joint_', tfr.sim.dir=NULL), 
						verbose=FALSE) {
	inp <- load.inputs(pop.pred$function.inputs, pop.pred$inputs$start.year, pop.pred$inputs$present.year, pop.pred$inputs$end.year, 
								pop.pred$wpp.year, verbose=verbose)
	for (item in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MIGm', 'MIGf', 'SRB', 'PASFR', 'MIGtype'))
		inp[[item]] <- NULL
	aggregated.countries <- list()
	aggr.obs.data <- list(male=NULL, female=NULL)
	if(verbose) cat('\nAggregating inputs using regional TFR and e0.')
	status.for.gui <- paste('out of', length(regions), 'regions.')
	gui.options <- list()
	obs.data <- inp$pop.matrix
	region.counter <- 0
	for(id in regions) {
		if(getOption('bDem.PopAgpred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			region.counter <- region.counter + 1
			gui.options$bDem.PopAgpred.status <- paste('finished', region.counter, status.for.gui)
			unblock.gtk('bDem.PopAgpred', gui.options)
		}

		countries <- get.countries.for.region(id, pop.pred)
		if(length(countries)==0) next
		inipop <- .aggregate.initial.pop(pop.pred, countries, id)
		for(item in c('POPm0', 'POPf0')) inp[[item]] <- rbind(inp[[item]], inipop[[item]])
		mort <- .aggregate.mortality.rates(pop.pred, countries, id)
		for(item in c('MXm', 'MXf')) inp[[item]] <- rbind(inp[[item]], mort[[item]])
		mig <- .aggregate.migration(pop.pred, countries, id)
		for(item in c('MIGm', 'MIGf')) inp[[item]] <- rbind(inp[[item]], mig[[item]])
		countries.index <- which(is.element(pop.pred$countries[,'code'], countries))
		inp$SRB <- rbind(inp$SRB, .aggregate.srb(pop.pred, countries, countries.index, id))
		inp$PASFR <- rbind(inp$PASFR, .aggregate.pasfr(pop.pred, countries, countries.index, id))
		inp$MIGtype <- rbind(inp$MIGtype, .aggregate.migtype(pop.pred, countries, countries.index, id))
		aggregated.countries[[as.character(id)]] <- countries
		# aggregate pop.matrix
		for(cidx in 1:length(countries.index)) {
			country.obs.idx <- grep(paste('^', countries[cidx], '_', sep=''), rownames(obs.data[['male']]), value=FALSE)
			if(cidx == 1) {
				aggr.obs.dataM <- obs.data[['male']][country.obs.idx,]
				aggr.obs.dataF <- obs.data[['female']][country.obs.idx,]
				rownames(aggr.obs.dataM) <- rownames(aggr.obs.dataF) <- sub(paste(countries[cidx], '_', sep=''), 
													paste(id, '_', sep=''), rownames(obs.data[['male']][country.obs.idx,]))
			} else {
				aggr.obs.dataM <- aggr.obs.dataM + obs.data[['male']][country.obs.idx,]
				aggr.obs.dataF <- aggr.obs.dataF + obs.data[['female']][country.obs.idx,]
			}
		}
		aggr.obs.data[['male']] <- rbind(aggr.obs.data[['male']], aggr.obs.dataM)
		aggr.obs.data[['female']] <- rbind(aggr.obs.data[['female']], aggr.obs.dataF)
	}
	dir <- if(is.null(inputs$tfr.sim.dir)) pop.pred$function.inputs$tfr.sim.dir else inputs$tfr.sim.dir
	if(!is.null(dir)) inp$TFRpred <- get.tfr.prediction(dir)
	dir <- if(is.null(inputs$e0F.sim.dir)) pop.pred$function.inputs$e0F.sim.dir else inputs$e0F.sim.dir
	if(!is.null(dir)) inp$e0Fpred <- get.e0.prediction(dir)
	if((is.null(inputs$e0M.sim.dir) || inputs$e0M.sim.dir == 'joint_') && has.e0.jmale.prediction(inp$e0Fpred)) {
		inp$e0Mpred <- get.e0.jmale.prediction(inp$e0Fpred)
	} else inp$e0Mpred <- get.e0.prediction(if(is.null(inputs$e0M.sim.dir)) pop.pred$function.inputs$e0M.sim.dir else inputs$e0M.sim.dir)
	outdir <- gsub('predictions', paste('aggregations', name, sep='_'), pop.output.directory(pop.pred))
	if(file.exists(outdir)) unlink(outdir, recursive=TRUE)
	dir.create(outdir, recursive=TRUE)
	aggr.pred <- do.pop.predict(regions, inp=inp, outdir=outdir, nr.traj=pop.pred$nr.traj, ages=pop.pred$ages, verbose=verbose)
	aggr.pred <- .cleanup.pop.before.save(aggr.pred, remove.cache=TRUE)
	aggr.pred$aggregation.method <- 'region'
	aggr.pred$aggregated.countries <- aggregated.countries
	aggr.pred$inputs$pop.matrix <- aggr.obs.data
	bayesPop.prediction <- aggr.pred
	save(bayesPop.prediction, file=file.path(outdir, 'prediction.rda'))
	return(bayesPop.prediction)
}

.aggregate.by.sum <- function(pop.pred, countries, what.names, id) {
	.sum.fcn <- function(age, index, what) {
		age.idx <- index & (pop.pred$inputs[[what]][,'age'] == age)
		return(colSums(pop.pred$inputs[[what]][age.idx,3:ncol.inp, drop=FALSE]))
	} 
	l <- 21	
	all.ages <- seq(0, by=5, length=l)
	age.labels <- c(paste(all.ages[1:(l-1)], '-', all.ages[2:l]-1, sep=''), paste(all.ages[l],'+',sep=''))
	res <- list()
	for(item in what.names) {
		ncol.inp <- ncol(pop.pred$inputs[[item]])
		idx <- is.element(pop.pred$inputs[[item]][,'country_code'], countries)
		values <- lapply(age.labels, .sum.fcn, index=idx, what=item)
		values <- matrix(unlist(values), nrow=length(values), byrow=TRUE)
		colnames(values) <- colnames(pop.pred$inputs[[item]])[3:ncol.inp]
		res[[item]] <- data.frame(country_code=rep(id,l), age=age.labels, values, check.names=FALSE)
	}
	return(res)
}

.aggregate.initial.pop <- function(pop.pred, countries, aggr.id) 
	return(.aggregate.by.sum(pop.pred, countries, c('POPm0', 'POPf0'), aggr.id))

.aggregate.migration <- function(pop.pred, countries, aggr.id) 
	return(.aggregate.by.sum(pop.pred, countries, c('MIGm', 'MIGf'), aggr.id))

.aggregate.mortality.rates <- function(pop.pred, countries, aggr.id) {
	popmatrix.all <- list(MXm=pop.pred$inputs$pop.matrix$male, MXf=pop.pred$inputs$pop.matrix$female)
	res <- list()
	for(item in c('MXm', 'MXf')) {
		mort.idx <- is.element(pop.pred$inputs[[item]][,'country_code'], countries)
		popmatrix <- popmatrix.all[[item]]
		pop.countries.ages <- matrix(unlist(strsplit(rownames(popmatrix), '_')), ncol=2, byrow=TRUE)
		pop.idx <- is.element(pop.countries.ages[,1], countries)
		years.mx <- unlist(strsplit(colnames(pop.pred$inputs[[item]])[3:ncol(pop.pred$inputs[[item]])], '-'))
		pop.colidx <- which(is.element(as.integer(colnames(popmatrix)), years.mx[seq(2,length(years.mx),by=2)]))
		res[[item]] <- as.data.frame(matrix(NA, nrow=22, ncol=ncol(pop.pred$inputs[[item]]), 
						dimnames=list(c(0, 1, seq(5, 95, by=5), '100+'), colnames(pop.pred$inputs[[item]]))))
		trim.age <- gsub(' ', '', pop.pred$inputs[[item]][,'age'])
		trim.age.unique <- unique(trim.age)
		for(age in rownames(res[[item]])) {
			mort.age <- if(age == '100+') '100' else age # in wpp2012 there is no '100+' (in the new one there is) [&& !(age %in%  trim.age.unique)]
			mort.age.idx <- mort.idx & trim.age == mort.age
			pop.age.idx <- rep(0, nrow(pop.countries.ages))
			if(age == "0" || age =="1") pattern <- '^0-4'
			else {
				if(age == "100+") pattern <- age
				else pattern <- paste('^', age, '-', sep='')
			}
			pop.age.idx[grep(pattern, pop.countries.ages[,2], value=FALSE)] <- TRUE
			pop.age.idx <- pop.idx & pop.age.idx
			#if(sum(pop.age.idx) == 0 || sum(mort.age.idx) == 0) next
			person.years <- popmatrix[pop.age.idx, pop.colidx] # doesn't need to be multiplied by years because it cancels out in the ratio below
			deaths <- colSums(pop.pred$inputs[[item]][mort.age.idx,3:ncol(pop.pred$inputs[[item]])] * person.years)
			res[[item]][age,3:ncol(pop.pred$inputs[[item]])] <- deaths/colSums(person.years)
		}
		res[[item]][,"age"] <- rownames(res[[item]])
		res[[item]][,"country_code"] <- aggr.id
	}
	return(res)
}

.aggregate.srb <- function(pop.pred, countries, countries.index, aggr.id) {
	pred.meansM <- pop.pred$quantilesMage[countries.index,1, "0.5", 2:dim(pop.pred$quantilesMage)[4]]
	pred.meansF <- pop.pred$quantilesFage[countries.index,1, "0.5", 2:dim(pop.pred$quantilesFage)[4]]
	aggr.srb <- data.frame(country_code=aggr.id, matrix(colSums(pred.meansM)/colSums(pred.meansF), nrow=1))
	colnames(aggr.srb) <- colnames(pop.pred$inputs[["SRB"]])
	return(aggr.srb)
}

.aggregate.pasfr <- function(pop.pred, countries, countries.index, aggr.id) {
	pred.meansF <- pop.pred$quantilesFage[countries.index,4:11, "0.5", 2:dim(pop.pred$quantilesFage)[4]]
	idx <- is.element(pop.pred$inputs[["PASFR"]][,'country_code'], countries)
	res.ncol <- ncol(pop.pred$inputs[["PASFR"]])
	res <- as.data.frame(matrix(NA, nrow=7, ncol=res.ncol, 
					dimnames=list(NULL, colnames(pop.pred$inputs[["PASFR"]]))))
	ages <- pop.pred$inputs[["PASFR"]][1:7,"age"]
	res[,1] <- aggr.id
	res[,2] <- ages
	for (iage in 1:7) {
		age.idx <- idx & pop.pred$inputs[["PASFR"]][,'age'] == ages[iage]
		age.pasfr <- pop.pred$inputs[["PASFR"]][age.idx,3:res.ncol]
		res[iage,3:res.ncol] <- colSums(age.pasfr * pred.meansF[,iage,])/colSums(pred.meansF[,iage,])
	}
	return(res)
}

.aggregate.migtype <- function(pop.pred, countries, countries.index, aggr.id) {
	idx <- is.element(pop.pred$inputs[["MIGtype"]][,'country_code'], countries)
	obsN <- round(pop.pred$quantiles[countries.index, "0.5", 1],0) # observed in present year
	res <- matrix(NA, nrow=1, ncol=3, dimnames=list(NULL, colnames(pop.pred$inputs[["MIGtype"]])))
	res[,1] <- aggr.id
	res[,2] <- max(pop.pred$inputs[["MIGtype"]][idx,"ProjFirstYear"])
	f <- factor(rep(pop.pred$inputs[["MIGtype"]][idx,"MigCode"], obsN), levels=c(0,9))
	res[,3] <- as.integer(f[which.max(tapply(rep(1,sum(obsN)), f, sum))])
	return(res)
}


pop.aggregate.countries <- function(pop.pred, regions, name, verbose=verbose, adjust=FALSE) {
	if(verbose) cat('\nAggregating using countries as inputs.')
	nreg <- length(regions)
	quantiles.to.keep <- as.numeric(dimnames(pop.pred$quantiles)[[2]])
	quant <- quantM <- quantF <- array(NA, c(nreg, dim(pop.pred$quantiles)[2:3]), dimnames=c(list(regions), dimnames(pop.pred$quantiles)[2:3]))
	quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(nreg, dim(pop.pred$quantilesMage)[2:4]),
						dimnames=c(list(regions), dimnames(pop.pred$quantilesMage)[2:4]))
	mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(nreg,dim(pop.pred$traj.mean.sd)[2:3]))
	obs.data <- pop.pred$inputs$pop.matrix
	aggr.obs.data <- list(male=NULL, female=NULL)
	outdir <- gsub('predictions', paste('aggregations', name, sep='_'), pop.output.directory(pop.pred))
	if(file.exists(outdir)) unlink(outdir, recursive=TRUE)
	dir.create(outdir, recursive=TRUE)
	aggregated.countries <- list()
	id.idx <- 0
	valid.regions <- rep(FALSE, length(regions))
	status.for.gui <- paste('out of', nreg, 'regions.')
	gui.options <- list()
	has.vital.events <- FALSE
	aggr.quantities <- c('totp', 'totpm', 'totpf', 'totp.hch', 'totpm.hch', 'totpf.hch')
	aggr.quantities.ve <- c('btm', 'btf', 'btm.hch', 'btf.hch', 
								'deathsm', 'deathsf', 'deathsm.hch', 'deathsf.hch',
								'migm', 'migf')
	# The next two lines need to be there for the R checker to get to know these objects
	totp <- totpm <- totpf <- totp.hch <- totpm.hch <- totpf.hch <- NULL
	btm <- btf <- deathsm <- deathsf <- migm <- migf <- btm.hch <- btf.hch <- deathsm.hch <- deathsf.hch <- NULL
	aggr.quantities.all <- aggr.quantities
	for(reg.idx in 1:length(regions)) {
		if(getOption('bDem.PopAgpred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			gui.options$bDem.PopAgpred.status <- paste('finished', reg.idx, status.for.gui)
			unblock.gtk('bDem.PopAgpred', gui.options)
		}
		id <- regions[reg.idx]
		if(verbose) cat('\nAggregating region ', id)
		countries <- get.countries.for.region(id, pop.pred)
		if(length(countries)==0) next
		id.idx <- id.idx + 1
		valid.regions[reg.idx] <- TRUE
		countries.index <- which(is.element(pop.pred$countries[,'code'], countries))
		e <- new.env()
		if(adjust && is.null(pop.pred$adjust.env)) pop.pred$adjust.env <- new.env()
		
		for(cidx in 1:length(countries.index)) {
			country.obs.idx <- grep(paste('^', countries[cidx], '_', sep=''), rownames(obs.data[['male']]), value=FALSE)
			traj.file <- file.path(pop.output.directory(pop.pred), paste('totpop_country', countries[cidx], '.rda', sep=''))
			load(traj.file, envir=e)
			if(adjust) adjust.trajectories(countries[cidx], e, pop.pred, pop.pred$adjust.env)
			if(has.vital.events || cidx == 1) {
				ve.file <- file.path(pop.output.directory(pop.pred), 
									paste('vital_events_country', countries[cidx], '.rda', sep=''))
				if(file.exists(ve.file)) {
					if(cidx == 1) {
						has.vital.events <- TRUE
						aggr.quantities.all <- c(aggr.quantities, aggr.quantities.ve)
						observed <- list()
					}
					load(ve.file, envir=e)
				}
			}
			if(cidx == 1) {
				for(par in aggr.quantities.all)
					assign(par, e[[par]])
				aggr.obs.dataM <- obs.data[['male']][country.obs.idx,, drop=FALSE]
				aggr.obs.dataF <- obs.data[['female']][country.obs.idx,, drop=FALSE]
				rownames(aggr.obs.dataM) <- rownames(aggr.obs.dataF) <- sub(paste(countries[cidx], '_', sep=''), 
																paste(id, '_', sep=''), rownames(obs.data[['male']][country.obs.idx,, drop=FALSE]))
				if(has.vital.events) {
					for(par in aggr.quantities.ve)
						if(!is.null(e$observed[[par]]))
							observed[[par]] <- e$observed[[par]]
				}
				trajectory.indices <- e$trajectory.indices											
			} else {
				for(par in aggr.quantities.all)
					assign(par, get(par) + e[[par]])
				aggr.obs.dataM <- aggr.obs.dataM + obs.data[['male']][country.obs.idx,]
				aggr.obs.dataF <- aggr.obs.dataF + obs.data[['female']][country.obs.idx,]
				if(has.vital.events) {
					for(par in names(observed))
						observed[[par]] <- observed[[par]] + e$observed[[par]]
				}
			}
		}
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch, trajectory.indices,
			 file = file.path(outdir, paste0('totpop_country', id, '.rda')))
		if(has.vital.events) {
			# asfert
			# Attach the second last observed year of female population in order to get last observed fertility
			tmp <- aggr.obs.dataF[4:10,as.character(pop.pred$proj.years.pop[1]-5),drop=FALSE]
			tmp <- abind(tmp, NULL, along=3)
			popfwprev <- abind(tmp[,,rep(1,dim(totpf)[3]), drop=FALSE], totpf[4:10,,,drop=FALSE], along=2)
			asfert <- 2*(btm + btf)/(popfwprev[,-dim(popfwprev)[2],,drop=FALSE]+
										popfwprev[,-1,,drop=FALSE])
			popfwprev.hch <- abind(tmp[,,rep(1,dim(totpf.hch)[3]), drop=FALSE], totpf.hch[4:10,,,drop=FALSE], along=2)
			asfert.hch <- 2*(btm.hch + btf.hch)/(popfwprev.hch[,-dim(popfwprev.hch)[2],,drop=FALSE]+
										popfwprev.hch[,-1,,drop=FALSE])
			# pasfert
			tfr <- apply(asfert, c(2,3), sum)
			pasfert <- asfert/abind(tfr, NULL, along=0)[rep(1,dim(asfert)[1]),,,drop=FALSE]*100
			tfr.hch <- apply(asfert.hch, c(2,3), sum)
			pasfert.hch <- asfert.hch/abind(tfr.hch, NULL, along=0)[rep(1,dim(asfert.hch)[1]),,,drop=FALSE]*100
			
			# TODO: mxm, mxf, mxm.hch, mxf.hch
			
			# asfert, pasfert for observed data
			observed <- within(observed, {
				tmp <- abind(aggr.obs.dataF[4:10, , drop=FALSE], NULL, along=3)
				if(dim(tmp)[2] > dim(btf)[2]+1) # if dimension of births doesn't match population
					tmp <- tmp[,-(1:(dim(tmp)[2]-dim(btf)[2]-1)),, drop=FALSE]
				asfert <- 2*(btm + btf)/(tmp[,-dim(tmp)[2],,drop=FALSE] + tmp[,-1,,drop=FALSE])
				tfr <- apply(asfert, c(2,3), sum)
				pasfert <- asfert/abind(tfr, NULL, along=0)[rep(1,dim(asfert)[1]),,,drop=FALSE]*100
				rm(tmp, tfr)
				# TODO: mxm, mxf, mxm.hch, mxf.hch
			})
			save(btm, btf, deathsm, deathsf, migm, migf, asfert, pasfert,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, pasfert.hch,
				observed, file=file.path(outdir, paste0('vital_events_country', id, '.rda')))
		}		
		quant[id.idx,,] = apply(totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sd[id.idx,1,] <- apply(totp, 1, mean, na.rm = TRUE)
		mean_sd[id.idx,2,] = apply(totp, 1, sd, na.rm = TRUE)
		
		for (i in 1:dim(totpm)[1]) {
			if(dim(totpm)[3] == 1) { # 1 trajectory
				quantMage[id.idx,i,,] <- matrix(totpm[i,,1], nrow=dim(quantMage)[3], ncol=dim(quantMage)[4], byrow=TRUE)
				quantFage[id.idx,i,,] <- matrix(totpf[i,,1], nrow=dim(quantFage)[3], ncol=dim(quantFage)[4], byrow=TRUE)
				quantPropMage[id.idx,i,,] <- matrix(totpm[i,,1]/totp[,1], nrow=dim(quantPropMage)[3], 
													ncol=dim(quantPropMage)[4], byrow=TRUE)
				
				quantPropFage[id.idx,i,,] <- matrix(totpf[i,,1]/totp[,1], nrow=dim(quantPropFage)[3], 
													ncol=dim(quantPropFage)[4], byrow=TRUE)
			} else { # multiple trajectories
				quantMage[id.idx,i,,] <- apply(totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantFage[id.idx,i,,] <- apply(totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropMage[id.idx,i,,] <- apply(totpm[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropFage[id.idx,i,,] <- apply(totpf[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			}
		}
		sstotpm <- colSums(totpm)
		quantM[id.idx,,] = apply(sstotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdM[id.idx,1,] <- apply(sstotpm, 1, mean, na.rm = TRUE)
		mean_sdM[id.idx,2,] = apply(sstotpm, 1, sd, na.rm = TRUE)
		sstotpf <- colSums(totpf)
		quantF[id.idx,,] = apply(sstotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdF[id.idx,1,] <- apply(sstotpf, 1, mean, na.rm = TRUE)
		mean_sdF[id.idx,2,] = apply(sstotpf, 1, sd, na.rm = TRUE)
		aggregated.countries[[as.character(id)]] <- countries
		aggr.obs.data[['male']] <- rbind(aggr.obs.data[['male']], aggr.obs.dataM)
		aggr.obs.data[['female']] <- rbind(aggr.obs.data[['female']], aggr.obs.dataF)
	}
	aggr.pred <- pop.pred
	which.reg.index <- function(x, set) return(which(set == x))
	reg.idx <- sapply(regions[valid.regions], which.reg.index, set=UNlocations[,'country_code']) 
	aggr.pred$countries=data.frame(code=UNlocations[reg.idx, 'country_code'], name=UNlocations[reg.idx, 'name'])
	aggr.pred$quantiles <- quant
	aggr.pred$quantilesM <- quantM
	aggr.pred$quantilesF <- quantF
	aggr.pred$quantilesMage <- quantMage
	aggr.pred$quantilesFage <- quantFage
	aggr.pred$quantilesPropMage <- quantPropMage
	aggr.pred$quantilesPropFage <- quantPropFage
	aggr.pred$traj.mean.sd <- mean_sd
	aggr.pred$traj.mean.sdM <- mean_sdM
	aggr.pred$traj.mean.sdF <- mean_sdF
	aggr.pred$aggregation.method <- 'country'
	aggr.pred$aggregated.countries <- aggregated.countries
	aggr.pred$inputs <- as.environment(as.list(pop.pred$inputs, all.names=TRUE)) # clone environment
	aggr.pred$inputs$pop.matrix <- aggr.obs.data
	bayesPop.prediction <- .cleanup.pop.before.save(aggr.pred, remove.cache=TRUE)
	save(bayesPop.prediction, file=file.path(outdir, 'prediction.rda'))
	cat('\nAggregations stored into', outdir, '\n')
	return(bayesPop.prediction)
}