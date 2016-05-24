if(getRversion() >= "2.15.1") utils::globalVariables(c("UNlocations", "MLTbx"))

pop.predict <- function(end.year=2100, start.year=1950, present.year=2015, wpp.year=2015,
						countries=NULL, output.dir = file.path(getwd(), "bayesPop.output"),
						inputs=list(
							popM=NULL,
							popF=NULL,
							mxM=NULL,
							mxF=NULL,
							srb=NULL,
							pasfr=NULL,
							mig.type=NULL,
							migM=NULL,
							migF=NULL,	
							e0F.file=NULL, e0M.file=NULL, 
							tfr.file=NULL, 
							e0F.sim.dir=NULL, e0M.sim.dir=NULL, 
							tfr.sim.dir=NULL,
							migMtraj=NULL, migFtraj=NULL	
						), nr.traj = 1000, keep.vital.events=FALSE,
						fixed.mx=FALSE, fixed.pasfr=FALSE, my.locations.file = NULL, 
						replace.output=FALSE, verbose=TRUE) {
	prediction.exist <- FALSE
	ages=seq(0, by=5, length=27)
	unblock.gtk.if.needed('reading inputs')
	if(!is.null(my.locations.file)) {
		UNlocations <- NULL # needed for R check not to complain
		UNlocations <<- read.delim(file=my.locations.file, comment.char='#', check.names=FALSE)
		if(verbose) cat('Loading ', my.locations.file, '.\n')
	} else bayesTFR:::load.bdem.dataset('UNlocations', wpp.year, envir=globalenv(), verbose=verbose)
	if(is.null(countries)) {
		if(!replace.output && has.pop.prediction(sim.dir=output.dir))
			stop('Prediction in ', output.dir,
				' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		inp <- load.inputs(inputs, start.year, present.year, end.year, wpp.year, fixed.mx=fixed.mx, verbose=verbose)
	}else {
		if(has.pop.prediction(output.dir) && !replace.output) {
			pred <- get.pop.prediction(output.dir)
			inp <- load.inputs(pred$function.inputs, pred$inputs$start.year, pred$inputs$present.year, pred$inputs$end.year, 
								pred$wpp.year, fixed.mx=pred$inputs$fixed.mx, all.countries=FALSE, 
								existing.mig=list(MIGm=pred$inputs$MIGm, MIGf=pred$inputs$MIGf, 
												obsMIGm=pred$inputs$observed$MIGm, obsMIGf=pred$inputs$observed$MIGf),
								verbose=verbose)
			if(!missing(inputs)) 
				warning('Projection already exists. Using inputs from existing projection. Use replace.output=TRUE for updating inputs.')
			nr.traj <- pred$nr.traj
			ages <- pred$ages
			prediction.exist <- TRUE
		} else inp <- load.inputs(inputs, start.year, present.year, end.year, wpp.year, fixed.mx=fixed.mx, all.countries=FALSE, verbose=verbose)
	}

	outdir <- file.path(output.dir, 'predictions')
	if(!prediction.exist) {
		if(!replace.output && has.pop.prediction(sim.dir=output.dir))
			stop('Prediction in ', outdir,
			' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		unlink(outdir, recursive=TRUE)
		.remove.cache.file(outdir)
	} else pop.cleanup.cache(pred)
	if(!is.null(countries) && is.na(countries[1])) { # all countries that are not included in the existing prediction
		all.countries <- intersect(unique(inp$POPm0[,'country_code']), UNcountries())
		country.codes <- if(!prediction.exist) all.countries
						else all.countries[!is.element(all.countries, pred$countries[,'code'])]
	} else {
		if(!is.null(countries)) {
			if (is.character(countries)) { # at least one of the codes is a character string
				for (icountry in 1:length(countries)) {
					if (is.character(countries[icountry])) {
						country.idx <- which(UNlocations[,'name'] == countries[icountry])
						if(length(country.idx) > 0)
							countries[icountry] <- UNlocations[country.idx,'country_code']
					}
				}
			}
			country.codes <- as.integer(countries)
		} else
			country.codes <- intersect(unique(inp$POPm0[,'country_code']), UNcountries())
	}
	do.pop.predict(country.codes, inp, outdir, nr.traj, ages, pred=if(prediction.exist) pred else NULL,
					keep.vital.events=keep.vital.events, fixed.mx=inp$fixed.mx, fixed.pasfr=fixed.pasfr, 
					function.inputs=inputs, verbose=verbose)
	invisible(get.pop.prediction(output.dir))
}

do.pop.predict <- function(country.codes, inp, outdir, nr.traj, ages, pred=NULL, keep.vital.events=FALSE, fixed.mx=FALSE, 
							fixed.pasfr=FALSE, function.inputs=NULL, verbose=FALSE) {
	not.valid.countries.idx <- c()
	countries.idx <- rep(NA, length(country.codes))

	for(icountry in 1:length(country.codes)) {
		country.idx <- which(UNlocations[,'country_code'] == country.codes[icountry])
		if(length(country.idx) == 0) {
			not.valid.countries.idx <- c(not.valid.countries.idx, icountry)
			next
		}
		countries.idx[icountry] <- country.idx
	}
	if(length(not.valid.countries.idx) > 0) {
		warning('Countries ', paste(country.codes[not.valid.countries.idx], collapse=', '), 
					' not found in the UNlocations dataset.')
		country.codes <- country.codes[-not.valid.countries.idx]
		countries.idx <- countries.idx[-not.valid.countries.idx]
	}
	ncountries <- length(country.codes)
	nr_project <- length(inp$proj.years)
	nages <- length(ages)
	if(!file.exists(outdir)) 
		dir.create(outdir, recursive=TRUE)
	present.and.proj.years <- c(inp$estim.years[length(inp$estim.years)], inp$proj.years)
	present.and.proj.years.pop <- present.and.proj.years + 2
	prediction.file <- file.path(outdir, 'prediction.rda')	
	quantiles.to.keep <- get.quantiles.to.keep()
	nquant <- length(quantiles.to.keep)
	PIs_cqp <- quantM <- quantF <- array(NA, c(ncountries, nquant, nr_project+1),
						dimnames=list(country.codes, quantiles.to.keep, present.and.proj.years.pop))
	quantMage <- quantFage <- quantPropMage <- quantPropFage <- array(NA, c(ncountries, nages, nquant, nr_project+1),
						dimnames=list(country.codes, ages, quantiles.to.keep, present.and.proj.years.pop))
	mean_sd <- mean_sdM <- mean_sdF <- array(NA, c(ncountries, 2, nr_project+1), 
						dimnames=list(country.codes, c('mean', 'sd'), present.and.proj.years.pop))
	mx.ages <- c(0,1,ages[2:nages])
	status.for.gui <- paste('out of', ncountries, 'countries.')
	gui.options <- list()
	inp.to.save <- list()
	# remove big or redundant items from inputs to be saved
	for(item in ls(inp)[!grepl('^migMpred$|^migFpred$|^TFRpred$|^e0Fpred$|^e0Mpred$|^estim.years$|^proj.years$|^wpp.years$', ls(inp))]) 
		inp.to.save[[item]] <- get(item, inp)
		
	for(cidx in 1:ncountries) {
		unblock.gtk.if.needed(paste('finished', cidx, status.for.gui), gui.options)
		country <- country.codes[cidx]
		country.idx <- countries.idx[cidx]
		if(verbose)
			cat('\nProgress: ', round((cidx-1)/ncountries * 100), '%; now processing ', country, ' ', as.character(UNlocations[country.idx,'name']))
		# Extract the country-specific stuff from the inputs
		inpc <- get.country.inputs(country, inp, nr.traj, UNlocations[country.idx,'name'])
		if(is.null(inpc)) next
		nr.traj <- min(ncol(inpc$TFRpred), nr.traj)		
		if(verbose)
			cat(' (', nr.traj, ' trajectories )')
		migr.modified <- .set.inp.migration.if.needed(inp, inpc, country)
		if(migr.modified) {
			for(what.mig in c('MIGm', 'MIGf')) {
				inp.to.save[[what.mig]] <- inp[[what.mig]]
				inp.to.save$observed[[what.mig]] <- inp$observed[[what.mig]]
			}
		}
		npred <- min(nrow(inpc$TFRpred), nr_project)
		npredplus1 <- npred+1
		totp <- matrix(NA, nrow=npredplus1, ncol=nr.traj, 
					dimnames=list(present.and.proj.years.pop, NULL))
		totpm <- totpf <- array(NA, dim=c(27, npredplus1, nr.traj), 
							dimnames=list(ages, present.and.proj.years.pop, NULL))
		nvariants <- nrow(inpc$TFRhalfchild)
		totp.hch <- matrix(NA, nrow=npredplus1, ncol=nvariants,
					dimnames=list(present.and.proj.years.pop, NULL))
		totpm.hch <- totpf.hch <- array(NA, dim=c(27, npredplus1, nvariants), 
					dimnames=list(ages, present.and.proj.years.pop, NULL))
		if(keep.vital.events) {
			btm <- btf <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm <- deathsf <- array(0, dim=c(27, npredplus1, nr.traj), dimnames=list(ages, present.and.proj.years, NULL))
			asfert <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			pasfert <- array(0, dim=c(7, npredplus1, nr.traj), dimnames=list(NULL, present.and.proj.years, NULL))
			btm.hch <- btf.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			deathsm.hch <- deathsf.hch <- array(0, dim=c(27, npredplus1, nvariants), dimnames=list(ages, present.and.proj.years, NULL))
			asfert.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			pasfert.hch <- array(0, dim=c(7, npredplus1, nvariants), dimnames=list(NULL, present.and.proj.years, NULL))
			mxm <- mxf <- array(0, dim=c(28, npredplus1, nr.traj), dimnames=list(mx.ages, present.and.proj.years, NULL))
			mxm.hch <- mxf.hch <- array(0, dim=c(28, npredplus1, nvariants), dimnames=list(mx.ages, present.and.proj.years, NULL))
			migMntraj <- if(is.null(inpc[['migMpred']])) 1 else dim(inpc[['migMpred']])[2]
			migFntraj <- if(is.null(inpc[['migFpred']])) 1 else dim(inpc[['migFpred']])[2]  
			migm <- array(0, dim=c(21, npredplus1, migMntraj), dimnames=list(ages[1:21], present.and.proj.years, NULL))
			migf <- array(0, dim=c(21, npredplus1, migFntraj), dimnames=list(ages[1:21], present.and.proj.years, NULL))
			# extend current mx to 27 age groups
			Kan.present <- KannistoAxBx.joint(inpc$MXm[,ncol(inpc$MXm), drop=FALSE], inpc$MXf[,ncol(inpc$MXf), drop=FALSE], compute.AxBx=FALSE)
		}
		debug <- FALSE
		#stop('')
		if(!fixed.mx) 
			MxKan <- runKannisto(inpc, inp$start.year, npred=npred) 
		else {
			MxKan <- runKannisto.noLC(inpc)
			LTres <- survival.fromLT(npred, MxKan, verbose=verbose, debug=debug)
		}
		#npasfr <- nrow(inpc$PASFR)
		if(keep.vital.events) observed <- compute.observedVE(inpc, inp$pop.matrix, inpc$MIGtype, MxKan, country, 
															estim.years=inp$estim.years)
		tfr.med <- apply(inpc$TFRpred, 1, median)[nrow(inpc$TFRpred)]
		for(itraj in 1:nr.traj) {
			if(any(is.na(inpc$TFRpred[,itraj]))) next
			if(!fixed.pasfr) 
				pasfr <- kantorova.pasfr(c(inpc$observed$TFRpred, inpc$TFRpred[,itraj]), inpc, 
										norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)
			else pasfr <- inpc$PASFR/100.
			asfr <- pasfr
			for(i in 1:nrow(pasfr)) asfr[i,] <- inpc$TFRpred[,itraj] * asfr[i,]
			if(!fixed.mx) LTres <- modifiedLC(npred, MxKan, inpc$e0Mpred[,itraj], 
									inpc$e0Fpred[,itraj], verbose=verbose, debug=debug)
			migpred <- list(M=NULL, F=NULL)
			for(sex in c('M', 'F')) {
				par <- paste0('mig', sex, 'pred')
				migpred[[sex]] <- as.matrix(if(is.null(inpc[[par]])) inpc[[paste0('MIG', tolower(sex))]] else inpc[[par]][,itraj,])
			}
			popres <- StoPopProj(npred, inpc, LTres, asfr, migpred, inpc$MIGtype, country.name=UNlocations[country.idx,'name'],
									keep.vital.events=keep.vital.events)
			totp[,itraj] <- popres$totpop
			totpm[,,itraj] <- popres$mpop
			totpf[,,itraj] <- popres$fpop
			if(keep.vital.events) {
				btm[,2:npredplus1,itraj] <- popres$mbt
				btm[1:dim(observed$btm)[1],1,itraj] <- observed$btm[,dim(observed$btm)[2],]
				btf[,2:npredplus1,itraj] <- popres$fbt
				btf[1:dim(observed$btf)[1],1,itraj] <- observed$btf[,dim(observed$btf)[2],]
				deathsm[,2:npredplus1,itraj] <- popres$mdeaths
				deathsm[1:dim(observed$deathsm)[1],1,itraj] <- observed$deathsm[,dim(observed$deathsm)[2],]
				deathsf[,2:npredplus1,itraj] <- popres$fdeaths
				deathsf[1:dim(observed$deathsf)[1],1,itraj] <- observed$deathsf[,dim(observed$deathsf)[2],]
				asfert[,2:npredplus1,itraj] <- asfr
				asfert[1:dim(observed$asfert)[1],1,itraj] <- observed$asfert[,dim(observed$asfert)[2],]
				pasfert[,2:npredplus1,itraj] <- pasfr*100
				pasfert[1:dim(pasfr)[1],1,itraj] <- inpc$observed$PASFR[,dim(inpc$observed$PASFR)[2]]				
				mxm[,2:npredplus1,itraj] <- LTres$mx[[1]]
				mxm[1:length(Kan.present$male$mx),1,itraj] <- Kan.present$male$mx
				mxf[,2:npredplus1,itraj] <- LTres$mx[[2]]
				mxf[1:length(Kan.present$female$mx),1,itraj] <- Kan.present$female$mx
				migtraj <- min(itraj, migMntraj)
				migm[,2:npredplus1,migtraj] <- migpred[['M']]
				migm[1:dim(inpc$observed$MIGm)[1],1,migtraj] <- inpc$observed$MIGm[,dim(inpc$observed$MIGm)[2]]
				migtraj <- min(itraj, migFntraj)
				migf[,2:npredplus1,migtraj] <- migpred[['F']]
				migf[1:dim(inpc$observed$MIGf)[1],1,migtraj] <- inpc$observed$MIGf[,dim(inpc$observed$MIGf)[2]]
			}
		}
		for (variant in 1:nvariants) { # compute the two half child variants
			if(!fixed.pasfr) 
				pasfr <- kantorova.pasfr(c(inpc$observed$TFRpred, inpc$TFRhalfchild[variant,]), inpc, 
										norms=inp$PASFRnorms, proj.years=inp$proj.years, 
										tfr.med=tfr.med)
			else pasfr <- inpc$PASFR/100.
			asfr <- pasfr
			for(i in 1:nrow(pasfr)) asfr[i,] <- inpc$TFRhalfchild[variant,] * asfr[i,]
			if(!fixed.mx) LTres <- modifiedLC(npred, MxKan, inpc$e0Mmedian, 
									inpc$e0Fmedian, verbose=verbose, debug=debug)
			migpred.hch <- list(M=NULL, F=NULL)
			for(sex in c('M', 'F')) {
				par <- paste0('mig', sex, 'median')
				migpred.hch[[sex]] <- as.matrix(if(is.null(inpc[[par]])) inpc[[paste0('MIG', tolower(sex))]] else inpc[[par]])
			}
			popres <- StoPopProj(npred, inpc, LTres, asfr, migpred.hch, inpc$MIGtype, 
							country.name=UNlocations[country.idx,'name'], 
							keep.vital.events=keep.vital.events)
			totp.hch[,variant] <- popres$totpop
			totpm.hch[,,variant] <- popres$mpop
			totpf.hch[,,variant] <- popres$fpop
			if(keep.vital.events) {
				btm.hch[,2:npredplus1,variant] <- popres$mbt
				btm.hch[,1,variant] <- btm[,1,1]
				btf.hch[,2:npredplus1,variant] <- popres$fbt
				btf.hch[,1,variant] <- btf[,1,1]
				deathsm.hch[,2:npredplus1,variant] <- popres$mdeaths
				deathsm.hch[,1,variant] <- deathsm[,1,1]
				deathsf.hch[,2:npredplus1,variant] <- popres$fdeaths
				deathsf.hch[,1,variant] <- deathsf[,1,1]
				asfert.hch[,2:npredplus1,variant] <- asfr
				asfert.hch[,1,variant] <- asfert[,1,1]
				pasfert.hch[,2:npredplus1,variant] <- pasfr*100
				pasfert.hch[,1,variant] <- pasfert[,1,1]
				mxm.hch[,2:npredplus1,variant] <- LTres$mx[[1]]
				mxf.hch[,2:npredplus1,variant] <- LTres$mx[[2]]
				mxm.hch[,1,variant] <- mxm[,1,1]
				mxf.hch[,1,variant] <- mxf[,1,1]
			}
		}
		trajectory.indices <- inpc$trajectory.indices
		save(totp, totpm, totpf, totp.hch, totpm.hch, totpf.hch, trajectory.indices,
			 file = file.path(outdir, paste0('totpop_country', country, '.rda')))
		if(keep.vital.events) 
			save(btm, btf, deathsm, deathsf, asfert, pasfert, mxm, mxf, migm, migf,
				btm.hch, btf.hch, deathsm.hch, deathsf.hch, asfert.hch, pasfert.hch, 
				mxm.hch, mxf.hch, 
				observed,
					file=file.path(outdir, paste0('vital_events_country', country, '.rda')))
		PIs_cqp[cidx,,] = apply(totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sd[cidx,1,] <- apply(totp, 1, mean, na.rm = TRUE)
		mean_sd[cidx,2,] = apply(totp, 1, sd, na.rm = TRUE)
		for (i in 1:nages) {
			if(nr.traj == 1) {
				quantMage[cidx,i,,] <- matrix(rep(totpm[i,,1],nquant) , nrow=nquant, byrow=TRUE)
				quantFage[cidx,i,,] <- matrix(rep(totpf[i,,1],nquant) , nrow=nquant, byrow=TRUE)
				quantPropMage[cidx,i,,] <- matrix(rep(totpm[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
				quantPropFage[cidx,i,,] <- matrix(rep(totpf[i,,1]/totp,nquant) , nrow=nquant, byrow=TRUE)
			} else {
				quantMage[cidx,i,,] <- apply(totpm[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantFage[cidx,i,,] <- apply(totpf[i,,], 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropMage[cidx,i,,] <- apply(totpm[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
				quantPropFage[cidx,i,,] <- apply(totpf[i,,]/totp, 1, quantile, quantiles.to.keep, na.rm = TRUE)
			}
		}
		stotpm <- colSums(totpm, na.rm=TRUE)
		quantM[cidx,,] = apply(stotpm, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdM[cidx,1,] <- apply(stotpm, 1, mean, na.rm = TRUE)
		mean_sdM[cidx,2,] = apply(stotpm, 1, sd, na.rm = TRUE)
		stotpf <- colSums(totpf, na.rm=TRUE)
		quantF[cidx,,] = apply(stotpf, 1, quantile, quantiles.to.keep, na.rm = TRUE)
		mean_sdF[cidx,1,] <- apply(stotpf, 1, mean, na.rm = TRUE)
		mean_sdF[cidx,2,] = apply(stotpf, 1, sd, na.rm = TRUE)
		
		#save updated meta file
		country.row <- UNlocations[country.idx,c('country_code', 'name')]
		colnames(country.row) <- c('code', 'name')
		

		
		if(!exists('bayesPop.prediction')) { # first pass
			bayesPop.prediction <- if(!is.null(pred)) .cleanup.pop.before.save(pred, remove.cache= country %in% pred$countries[,'code']) 
					else structure(list(
							nr.traj = nr.traj,	
							# assign empty arrays
							quantiles = PIs_cqp[c(),,,drop=FALSE],
               				traj.mean.sd = mean_sd[c(),,,drop=FALSE],
               				quantilesM = quantM[c(),,,drop=FALSE], 
               				traj.mean.sdM = mean_sdM[c(),,,drop=FALSE],
               				quantilesF = quantF[c(),,,drop=FALSE], 
               				traj.mean.sdF = mean_sdF[c(),,,drop=FALSE],
               				quantilesMage = quantMage[c(),,,,drop=FALSE], 
               				quantilesFage = quantFage[c(),,,,drop=FALSE], 
               				quantilesPropMage = quantPropMage[c(),,,,drop=FALSE], 
               				quantilesPropFage = quantPropFage[c(),,,,drop=FALSE],
               				estim.years=inp$estim.years, 
               				proj.years=present.and.proj.years, # includes present period (middle of periods)
               				proj.years.pop=present.and.proj.years.pop, # end of periods
               				wpp.year = inp$wpp.year,
			   				inputs = inp.to.save, # save as list because environment takes much more space
			   				function.inputs=function.inputs,
			   				countries=as.data.frame(matrix(NA, nrow=0, ncol=2, dimnames=list(NULL, c('code', 'name')))),
			   				ages=ages), class='bayesPop.prediction')
		}
		idx.in.pred.overwrite <- which(bayesPop.prediction$countries[,'code'] == country)
		if(length(idx.in.pred.overwrite)>0) {
			bayesPop.prediction$quantiles[idx.in.pred.overwrite,,] <- PIs_cqp[cidx,,,drop=FALSE]
			bayesPop.prediction$traj.mean.sd[idx.in.pred.overwrite,,] <- mean_sd[cidx,,,drop=FALSE]
			bayesPop.prediction$traj.mean.sdM[idx.in.pred.overwrite,,] <- mean_sdM[cidx,,,drop=FALSE]
			bayesPop.prediction$traj.mean.sdF[idx.in.pred.overwrite,,] <- mean_sdF[cidx,,,drop=FALSE]
			bayesPop.prediction$quantilesM[idx.in.pred.overwrite,,] <- quantM[cidx,,,drop=FALSE]
			bayesPop.prediction$quantilesF[idx.in.pred.overwrite,,] <- quantF[cidx,,,drop=FALSE]
			bayesPop.prediction$quantilesMage[idx.in.pred.overwrite,,,] <- quantMage[cidx,,,,drop=FALSE]
			bayesPop.prediction$quantilesFage[idx.in.pred.overwrite,,,] <- quantFage[cidx,,,,drop=FALSE]
			bayesPop.prediction$quantilesPropMage[idx.in.pred.overwrite,,,] <- quantPropMage[cidx,,,,drop=FALSE]
			bayesPop.prediction$quantilesPropFage[idx.in.pred.overwrite,,,] <- quantPropFage[cidx,,,,drop=FALSE]
		} else { 
			bayesPop.prediction$quantiles <- abind(bayesPop.prediction$quantiles, 
												PIs_cqp[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$traj.mean.sd <- abind(bayesPop.prediction$traj.mean.sd, 
												mean_sd[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$traj.mean.sdM <- abind(bayesPop.prediction$traj.mean.sdM, 
												mean_sdM[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$traj.mean.sdF <- abind(bayesPop.prediction$traj.mean.sdF, 
												mean_sdF[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesM <- abind(bayesPop.prediction$quantilesM, 
												quantM[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesF <- abind(bayesPop.prediction$quantilesF, 
												quantF[cidx,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesMage <- abind(bayesPop.prediction$quantilesMage, 
												quantMage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesFage <- abind(bayesPop.prediction$quantilesFage, 
												quantFage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesPropMage <- abind(bayesPop.prediction$quantilesPropMage, 
												quantPropMage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$quantilesPropFage <- abind(bayesPop.prediction$quantilesPropFage, 
												quantPropFage[cidx,,,,drop=FALSE], along=1)
			bayesPop.prediction$countries <- rbind(bayesPop.prediction$countries, country.row)
		}
		if(migr.modified)
			bayesPop.prediction$inputs <- inp.to.save
		save(bayesPop.prediction, file=prediction.file)
	} 
	cat('\nPrediction stored into', outdir, '\n')
	return(bayesPop.prediction)
}

read.pop.file <- function(file) 
	return(read.delim(file=file, comment.char='#', check.names=FALSE))
	
load.wpp.dataset <- function(...)
	bayesTFR:::load.bdem.dataset(...)
	
read.bayesPop.file <- function(file)
	return(get(do.call('data', list(strsplit(file, '.', fixed=TRUE)[[1]][-2]))))

load.inputs <- function(inputs, start.year, present.year, end.year, wpp.year, fixed.mx=FALSE, 
						all.countries=TRUE, existing.mig=NULL, verbose=FALSE) {
	observed <- list()
	pop.ini.matrix <- pop.ini <- list(M=NULL, F=NULL)
	# Get initial population counts
	for(sex in c('M', 'F')) {
		dataset.name <- paste0('pop', sex)
		if(is.null(inputs[[dataset.name]])) 
			POP0 <- load.wpp.dataset(dataset.name, wpp.year)	
		else POP0 <- read.pop.file(inputs[[dataset.name]])
		num.columns <- grep('^[0-9]{4}$', colnames(POP0), value=TRUE) # values of year-columns
		if(!is.element(as.character(present.year), colnames(POP0))) {
			stop('Wrong present.year. ', present.year, ' not available in the ', dataset.name, ' file.\nAvailable years: ',
				paste(num.columns, collapse=', '))
		}
		num.columns <- num.columns[which(as.integer(num.columns)<= present.year)]
		pop.ini.matrix[[sex]] <- POP0[,num.columns, drop=FALSE]
		dimnames(pop.ini.matrix[[sex]]) <- list(paste(POP0[,'country_code'], POP0[,'age'], sep='_'), 
									as.character(as.integer(num.columns)))
		pop.ini[[sex]] <- POP0[,c('country_code', 'age', present.year)]
	}
	POPm0 <- pop.ini[['M']]
	POPf0 <- pop.ini[['F']]

	# Get death rates
	MXm.pred <- MXf.pred <- NULL
	if(is.null(inputs$mxM)) 
		MXm <- load.wpp.dataset('mxM', wpp.year)
	else MXm <- read.pop.file(inputs$mxM)
	names.MXm.data <- names(MXm)
	num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.MXm.data) # index of year-columns
	cols.starty <- as.integer(substr(names.MXm.data[num.columns], 1,4))
    cols.endy <- as.integer(substr(names.MXm.data[num.columns], 6,9))
	start.index <- which((cols.starty <= start.year) & (cols.endy > start.year))
	present.index <- which((cols.endy >= present.year) & (cols.starty < present.year))
	#estim.periods <- names.MXm.data[num.columns[start.index:present.index]]
	estim.periods <- names.MXm.data[num.columns[1:present.index]]
	start.year <- cols.starty[start.index]
	
	if(fixed.mx) {
		end.index <- which((cols.endy >= end.year) & (cols.starty < end.year))
		proj.periods <- names.MXm.data[num.columns[(present.index+1):end.index]]
		MXm.pred <- MXm[,c('country_code', 'age', proj.periods)]
	}
	MXm <- MXm[,c('country_code', 'age', estim.periods)]
	if(is.null(inputs$mxF)) 
		MXf <- load.wpp.dataset('mxF', wpp.year)
	else MXf <- read.pop.file(inputs$mxF)
	if(fixed.mx) MXf.pred <- MXf[,c('country_code', 'age', proj.periods)]
	MXf <- MXf[,c('country_code', 'age', estim.periods)]
	
	estim.years <- cols.starty[start.index:present.index] + 3
	# Get sex ratio at birth
	if(is.null(inputs$srb)) 
		SRB <- load.wpp.dataset('sexRatio', wpp.year)
	else SRB <- read.pop.file(inputs$srb)
	names.SRB.data <- names(SRB)
	num.columns <- grep('^[0-9]{4}.[0-9]{4}$', names.SRB.data) # index of year-columns
	cols.starty <- as.integer(substr(names.SRB.data[num.columns], 1,4))
    cols.endy <- as.integer(substr(names.SRB.data[num.columns], 6,9))
	start.index <- which((cols.starty <= present.year) & (cols.endy > present.year))
	end.index <- which((cols.endy >= end.year) & (cols.starty < end.year))
	if(length(end.index) == 0) {
		end.index <- length(num.columns)
		warning('Data for SexRatioAtBirth not available for all projection periods.\nLast projection period set to ', 
					names.SRB.data[num.columns[end.index]])
	}
	proj.periods <- names.SRB.data[num.columns[start.index:end.index]]
	obs.periods <- NULL
	if(start.index > 1) {
		obs.periods <- names.SRB.data[num.columns[1:(start.index-1)]]
		observed$SRB <- SRB[,c('country_code', obs.periods)]
	}
	SRB <- SRB[,c('country_code', proj.periods)]
	proj.years <- cols.starty[start.index:end.index] + 3
	
	# Get percentage age-specific fertility rate
	if(is.null(inputs$pasfr)) 
		PASFR <- load.wpp.dataset('percentASFR', wpp.year)
	else PASFR <- read.pop.file(inputs$pasfr)
	if(!is.null(obs.periods)) {
		avail.obs.periods <- is.element(obs.periods, colnames(PASFR))
		observed$PASFR <- PASFR[,c('country_code', 'age', obs.periods[avail.obs.periods])]
	}
	PASFR <- PASFR[,c('country_code', 'age', proj.periods)]
	
	# Get migration type and base year
	if(is.null(inputs$mig.type)) 
		vwBase <- read.bayesPop.file(paste('vwBaseYear', wpp.year, '.txt', sep=''))
	else vwBase <- read.pop.file(inputs$mig.type)
	MIGtype <- vwBase[,c('country_code', 'ProjFirstYear', 'MigCode')]

	create.pattern <- function(dataset, columns) {
		pattern <- data.frame(dataset[,'country_code'])
		for(col in columns)
			if(col %in% colnames(dataset))
				pattern <- cbind(pattern, dataset[,col])
		if(ncol(pattern)==1) pattern <- NULL
		else colnames(pattern) <- c('country_code', columns)[1:ncol(pattern)]
		return(pattern)
	}
	MXpattern <- create.pattern(vwBase, c("AgeMortalityType", "AgeMortalityPattern", "LatestAgeMortalityPattern", "SmoothLatestAgeMortalityPattern", "WPPAIDS"))
	PASFRpattern <- create.pattern(vwBase, c("PasfrNorm", paste0("Pasfr", .remove.all.spaces(levels(vwBase$PasfrNorm)))))
	# Get age-specific migration
	wppds <- data(package=paste0('wpp', wpp.year))
	recon.mig <- NULL
	if(is.null(inputs[['migM']])) {# cannot be inputs$migM because it would take the value of inputs$migMpred
	 	if('migrationM' %in% wppds$results[,'Item'])
			MIGm <- load.wpp.dataset('migrationM', wpp.year)
		else {
			if(all.countries) { # reconstruct migration for all countries
				recon.mig <- age.specific.migration(wpp.year=wpp.year, verbose=verbose)
				MIGm <- recon.mig$male
			} else {
				if (!is.null(existing.mig)) { # Migration dataset already exists, e.g. from a priovous simulation for different country
					MIGm <- existing.mig$MIGm
				} else {
					# Projection for a subset of countries, therefore migration will be recontructed later (in get.country.inputs).
					# Here create only a dataframe filled with NAs (its structure is taken from the wpp2012 dataset) 
					MIGm <- load.wpp.dataset('migrationM', 2012)
					MIGm[,!colnames(MIGm) %in% c("country", "country_code", "age")] <- NA
				}
			}
		}
	} else MIGm <- read.pop.file(inputs[['migM']])
	if(is.null(inputs[['migF']])) {
		if('migrationF' %in% wppds$results[,'Item'])
			MIGf <- load.wpp.dataset('migrationF', wpp.year)
		else {
			if(all.countries) { # reconstruct migration for all countries
				if(is.null(recon.mig)) recon.mig <- age.specific.migration(wpp.year=wpp.year, verbose=verbose)
				else MIGf <- recon.mig$female
			} else { # see the comment for MIGm above
				if (!is.null(existing.mig)) { 
					MIGf <- existing.mig$MIGf
				} else {
					MIGf <- load.wpp.dataset('migrationF', 2012)
					MIGf[,!colnames(MIGm) %in% c("country", "country_code", "age")] <- NA
				}
			}
		}
	} else MIGf <- read.pop.file(inputs[['migF']])
	if(!is.null(obs.periods)) {
		if (!is.null(existing.mig)) { # Migration dataset already exists, e.g. from a priovous simulation for different country
			observed$MIGm <- existing.mig$obsMIGm
			observed$MIGf <- existing.mig$obsMIGf
		} else {
			avail.obs.periods <- is.element(obs.periods, colnames(MIGm))
			observed$MIGm <- MIGm[,c('country_code', 'age', obs.periods[avail.obs.periods])]
			observed$MIGf <- MIGf[,c('country_code', 'age', obs.periods[avail.obs.periods])]
		}
	}
	MIGm <- MIGm[,c('country_code', 'age', proj.periods)]
	MIGf <- MIGf[,c('country_code', 'age', proj.periods)]
	# Get migration trajectories if available
	migMpred <- migFpred <- NULL
	for(sex in c('M', 'F')) {
		if(is.null(inputs[[paste0('mig',sex,'traj')]])) next
		file.name <- inputs[[paste0('mig',sex,'traj')]]
		if(!file.exists(file.name))
			stop('File ', file.name, ' does not exist.')
		# comma separated trajectories file
		var.name <- paste0('mig',sex, 'pred')
		if(verbose) cat('\nLoading ', file.name)
		migpred.raw <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
		migpred <- migpred.raw[,c('LocID', 'Year', 'Trajectory', 'Age', 'Migration')]
		colnames(migpred) <- c('country_code', 'year', 'trajectory', 'age', 'value')
		assign(var.name, migpred)
	}
	
	# Get life expectancy
	e0F.wpp.median.loaded <- FALSE
	e0Fpred <- e0Mpred <- NULL
	if(!fixed.mx){
	if(!is.null(inputs$e0F.file)) { # female
		if(inputs$e0F.file == 'median_') {
			e0Fpred <- .load.wpp.traj('e0F', wpp.year, median.only=TRUE)
			e0F.wpp.median.loaded <- TRUE
		} else {
			file.name <-  inputs$e0F.file
			if(!file.exists(file.name))
				stop('File ', file.name, 
					' does not exist.\nSet e0F.sim.dir, e0F.file or change WPP year.')
		 	# comma separated trajectories file
		 	if(verbose) cat('\nLoading ', file.name)
			e0Fpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
			e0Fpred <- e0Fpred[,c('LocID', 'Year', 'Trajectory', 'e0')]
			colnames(e0Fpred) <- c('country_code', 'year', 'trajectory', 'value')
		}
	} else {
		if(!is.null(inputs$e0F.sim.dir)) { 
			if(inputs$e0F.sim.dir == 'median_') {
				e0Fpred <- .load.wpp.traj('e0F', wpp.year, median.only=TRUE)
				e0F.wpp.median.loaded <- TRUE
			} else 
				e0Fpred <- get.e0.prediction(inputs$e0F.sim.dir, mcmc.dir=NA)
		} else e0Fpred <- .load.wpp.traj('e0F', wpp.year)			
	}
	
	if(!is.null(inputs$e0M.file)) { # male
		if(inputs$e0M.file == 'median_')
			e0Mpred <- .load.wpp.traj('e0M', wpp.year, median.only=TRUE)
		else {
			file.name <-  inputs$e0M.file
			if(!file.exists(file.name)) 
				stop('File ', file.name, 
					' does not exist.\nSet e0M.sim.dir, e0M.file or change WPP year.')
			if(verbose) cat('\nLoading ', file.name)
			e0Mpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
			e0Mpred <- e0Mpred[,c('LocID', 'Year', 'Trajectory', 'e0')]
			colnames(e0Mpred) <- c('country_code', 'year', 'trajectory', 'value')
		}
	} else {
		if(!is.null(inputs$e0M.sim.dir)) { 
			if(inputs$e0M.sim.dir == 'joint_') {
				if(e0F.wpp.median.loaded) e0Mpred <- .load.wpp.traj('e0M', wpp.year)
				else {
					if(!has.e0.jmale.prediction(e0Fpred))
						stop('No joint prediction for female and male available. Correct the e0M.sim.dir argument.' )
					e0Mpred <- get.e0.jmale.prediction(e0Fpred)
				}
			} else e0Mpred <- get.e0.prediction(inputs$e0M.sim.dir, mcmc.dir=NA) # independent from female
		} else
			e0Mpred <- .load.wpp.traj('e0M', wpp.year)
	}
	}
	# Get TFR
	if(!is.null(inputs$tfr.file)) {
		if(inputs$tfr.file == 'median_')
			TFRpred <- .load.wpp.traj('tfr', wpp.year, median.only=TRUE)
		else {
			file.name <- inputs$tfr.file
			if(!file.exists(file.name))
				stop('File ', file.name, 
					' does not exist.\nSet tfr.sim.dir, tfr.file or change WPP year.')
			if(verbose) cat('\nLoading ', file.name, '\n')
			TFRpred <- read.csv(file=file.name, comment.char='#', check.names=FALSE)
			TFRpred <- TFRpred[,c('LocID', 'Year', 'Trajectory', 'TF')]
			colnames(TFRpred) <- c('country_code', 'year', 'trajectory', 'value')
		} 
	} else {
		if(!is.null(inputs$tfr.sim.dir)) 
			TFRpred <- get.tfr.prediction(inputs$tfr.sim.dir, mcmc.dir=NA)
		else TFRpred <- .load.wpp.traj('tfr', wpp.year)
	}
	
	inp <- new.env()
	for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MXm.pred', 'MXf.pred', 'MXpattern', 'SRB',
				'PASFR', 'PASFRpattern', 'MIGtype', 'MIGm', 'MIGf',
				'e0Mpred', 'e0Fpred', 'TFRpred', 'migMpred', 'migFpred', 'estim.years', 'proj.years', 'wpp.year', 
				'start.year', 'present.year', 'end.year', 'fixed.mx', 'observed'))
		assign(par, get(par), envir=inp)
	inp$pop.matrix <- list(male=pop.ini.matrix[['M']], female=pop.ini.matrix[['F']])
	inp$PASFRnorms <- compute.pasfr.global.norms(inp)
	return(inp)
}

.set.inp.migration.if.needed <- function(inputs, inpc, country) {
	# If the age-specific migration values in "inputs" are NA 
	# (e.g. because it was reconstructed later in the code), 
	# replace those with values from the country-specific inputs.
	migr.modified <- FALSE
	for(what.mig in c('MIGm', 'MIGf')) {
		cidx <- which(inputs[[what.mig]]$country_code==country)
		cols <- intersect(colnames(inputs[[what.mig]]), colnames(inpc[[what.mig]]))
		if(any(!is.na(inputs[[what.mig]][cidx,cols]))) next
		inputs[[what.mig]][cidx,cols] <- inpc[[what.mig]][,cols]
		cidx <- which(inputs$observed[[what.mig]]$country_code==country)
		cols <- intersect(colnames(inputs$observed[[what.mig]]), colnames(inpc$observed[[what.mig]]))
		inputs$observed[[what.mig]][cidx,cols] <- inpc$observed[[what.mig]][,cols]
		migr.modified <- TRUE
	}
	return(migr.modified)
}

.load.wpp.traj <- function(type, wpp.year, median.only=FALSE) {
	dataset.obs <- dataset.low <- dataset.high <- NA
	if(type %in% c('e0F', 'e0M')){
		if(wpp.year < 2010) stop('e0 projections not available for wpp 2008.')
		#if(wpp.year == 2010) 
		dataset <- paste(type, 'proj', sep='')
		#else { # wpp 2012
		#	dataset <- paste(type, 'projMed', sep='')
		#}
		dataset.obs <- type
	} else { # tfr
		if(wpp.year < 2010) dataset <- type
		else {
			dataset <- paste(type, 'projMed', sep='')
			dataset.obs <- type
			if(!median.only) {
				dataset.low <- paste(type, 'projLow', sep='')
				dataset.high <- paste(type, 'projHigh', sep='')
			}
		}
	}
	if(!is.na(dataset.obs)) 
		pred.obs <- bayesTFR:::load.bdem.dataset(dataset.obs, wpp.year)
		
	pred.all <- NULL
	itraj <- 1
	for(dataset.name in c(dataset, dataset.low, dataset.high)) {
		if(is.na(dataset.name)) next
		pred <- bayesTFR:::load.bdem.dataset(dataset.name, wpp.year)
		remove.cols <- which(colnames(pred) %in% c('name', 'country', 'last.observed'))
		pred <- pred[,-remove.cols]
		if(!is.na(dataset.obs)) {
			remove.cols <- which(colnames(pred.obs) %in% c('name', 'country', 'last.observed', 
								colnames(pred)[-which(colnames(pred)=='country_code')]))
			pred <- merge(pred.obs[,-remove.cols], pred, by='country_code')
		}
		ncols <- ncol(pred)
		nonnum.idx <- which(colnames(pred)=='country_code')
		cnames <- colnames(pred)[-nonnum.idx]
		colnames(pred)[-nonnum.idx] <- paste(type, cnames, sep='')
		pred.long <- reshape(pred, direction='long', varying=(1:ncols)[-nonnum.idx], v.names=type, times=cnames)
		pred.long <- cbind(pred.long, year=as.integer(substr(pred.long$time,1,4))+3, trajectory=itraj)
		pred.long <- pred.long[,c('country_code', 'year', 'trajectory', type)]
		pred.all <- rbind(pred.all, pred.long)
		itraj <- itraj + 1
	}
	colnames(pred.all) <- c('country_code', 'year', 'trajectory', 'value')
	return(pred.all)
}

.get.migration.traj <- function(pred, par, country) {
		cidx <- pred$inputs[[par]][,'country_code'] == country 
		idx <- cidx & is.element(pred$inputs[[par]][,'year'], pred$inputs$proj.years)
		if(sum(idx) == 0) return(NULL)
		migdf <- pred$inputs[[par]][idx,-1]
		utrajs <- sort(unique(migdf$trajectory))
		ntrajs <- length(utrajs)
		migdf$age <- gsub("^\\s+|\\s+$", "", migdf$age) # trim leading and trailing whitespace
		lyears <- length(pred$inputs$proj.years)
		sorted.df <- data.frame(year=rep(pred$inputs$proj.years, each=ntrajs*21), trajectory=rep(rep(utrajs, each=21), times=lyears),
									age=c(paste(seq(0,95,by=5), seq(4,99,by=5), sep='-'), '100+'))
		# this is to get rows of the data frame in a particular order
		migdf <- merge(sorted.df, migdf, sort=FALSE)
		res <- array(migdf$value, dim=c(21, ntrajs, lyears))
		dimnames(res) <- list(1:21,  NULL, pred$inputs$proj.years)
		return(res)	
}

.pasfr.norm.name <- function(norms)
	return(paste0('Pasfr', .remove.all.spaces(norms)))
	
compute.pasfr.global.norms <- function(inputs) {
	pattern <- inputs$PASFRpattern
	if(is.null(pattern) || !('PasfrNorm' %in% colnames(pattern))) return(NULL)
	norms <- .pasfr.norm.name(levels(pattern$PasfrNorm))
	result <- list()
	for(norm in norms) {
		tpasfr <- NULL
		countries <- pattern$country_code[which(pattern[[norm]]==1)]
		for(country in countries) {
			pasfr <- .get.par.from.inputs('PASFR', inputs$observed, country)
			tpasfr <- if(is.null(tpasfr)) pasfr else tpasfr + pasfr
		}
		tpasfr <- tpasfr/(length(countries)*100)
		result[[norm]] <- scale(tpasfr, center=FALSE, scale=colSums(tpasfr))*100
	}
	return(result)
}

kantorova.pasfr <- function(tfr, inputs, norms, proj.years, tfr.med) {
	logit <- function(x) log(x/(1-x))
	inv.logit <- function(x) exp(x)/(1+exp(x))
	compute.mac <- function(x) {
		factors <- seq(17.5, by=5, length=dim(x)[1])
		mac <- rep(0, dim(x)[2])
        for(iage in 1:dim(x)[1]) 
        	mac <- mac + x[iage,]*factors[iage]/100.
        return(mac)
	}
	update.by.mac <- function(x, sp3i) {
		if(sp3i >= dim(x)[2]) return(x)
		phase3i <- seq(sp3i, dim(x)[2])
		mac <- compute.mac(x[,phase3i]*100)		
		mac.norm <- compute.mac(matrix(gnorm, ncol=1))
		maxi <- which.max(mac)
		if(mac[maxi] <= mac.norm)
			return(x)
		x[,phase3i[maxi:length(phase3i)]] <- x[,phase3i[maxi]]
		#stop('')
		return(x)
	}
	pattern <- inputs$PASFRpattern
	if(is.null(pattern)) return(inputs$PASFR)
	min.value <- 1e-3	
	pasfr.obs <- inputs$observed$PASFR
	
	years <- as.integer(names(tfr))
	if(length(years)==0)
		years <- sort(seq(proj.years[length(proj.years)], length=length(tfr), by=-5))
	lyears <- length(years)
	years.long <- c(years, seq(years[lyears]+5, by=5, length=15))
	tobs <- lyears - length(proj.years)
	end.year <- years[lyears]
	end.phase2 <- bayesTFR:::find.lambda.for.one.country(tfr, lyears)
	start.phase3 <- end.phase2 + 1
	#stop('')
	if(start.phase3 > lyears) { # Phase 3 hasn't start (Case 2)
		if(tfr[lyears] > 1.8) { # regress the last four points to approximate start of Phase 3
			df <- data.frame(tfr=tfr[(lyears-3):lyears], time=(lyears-3):lyears)
			reg <- lm(tfr~time, df)
			if(reg$coefficients[2] < -1e-3) {# use only if it has negative slope 
				start.phase3 <- min(round((1.8-reg$coefficients[1])/reg$coefficients[2],0)+1, lyears+10)
			} else {
				start.phase3 <- lyears+10
			}			
		}
		#endT <- years.long[min(start.phase3+5, length(years.long))]
		endT <- years.long[start.phase3+5]
	} else { # Case 1
		smaller.than.median <- tfr[start.phase3:lyears] < tfr.med
		if(all(smaller.than.median)) { # t_u does not exist
			#endT <- years.long[max(start.phase3+10, tobs+10)]
			endT <- years.long[max(lyears, start.phase3+5)]
		} else { # t_u exists
			first.larger <- which(!smaller.than.median)[1] + start.phase3 - 1
			endT <- years.long[max(first.larger, tobs+2)]
		}
	}
	#endT <- years.long[max(start.phase3+5, tobs+5)] # no upper bound
	startTi <- which(years == proj.years[1])
	gnorm <- norms[[.pasfr.norm.name(pattern[,'PasfrNorm'])]]
	gnorm <- gnorm[, ncol(gnorm)] # global norm from the last time period 
	asfr1 <- asfr2 <- res.asfr <- matrix(0, nrow=length(gnorm), ncol=length(proj.years))
	t.r <- years[startTi-1]
	tau.denominator <- endT - t.r
	p.r <- pasfr.obs[,ncol(pasfr.obs)]/100. # last observed pasfr
	p.r <- pmax(p.r, min.value)
	p.r <- p.r/sum(p.r)
	logit.pr <- logit(p.r)
	logit.dif <- logit(gnorm/100.) - logit.pr
	for(t in 1:ncol(asfr1)){
		asfr1[,t] <- logit.pr + min((years[t+tobs] - t.r)/tau.denominator, 1)*logit.dif
	}
	asfr1 <- inv.logit(asfr1)
	asfr1 <-  scale(asfr1, center=FALSE, scale=colSums(asfr1))
	
	p.e <- pasfr.obs[,ncol(pasfr.obs)-2]/100.
	p.e <- pmax(p.e, min.value)
	p.e <- p.e/sum(p.e)
	if(startTi < 4) { # not enough observed data
		yd <- years[1] - 5 * (3-startTi)
	} else yd <- years[startTi-3]
	tau.denominator2 <- t.r - yd
	logit.dif <- logit.pr - logit(p.e)
	for(t in 1:ncol(asfr2)){
		asfr2[,t] <- logit.pr + ((years[t+tobs] - t.r)/tau.denominator2) *logit.dif
	}
	#stop('')
	asfr2 <- inv.logit(asfr2)
	asfr2 <-  scale(asfr2, center=FALSE, scale=colSums(asfr2))
	
	logit.asfr1 <- logit(asfr1)
	logit.asfr2 <- logit(asfr2)
	for(t in 1:ncol(res.asfr)){
		tau <- min((years[t+tobs] - t.r)/tau.denominator, 1)
		res.asfr[,t] <- tau*logit.asfr1[,t] + (1-tau)*logit.asfr2[,t]
	}
	res.asfr <- inv.logit(res.asfr)
	res.asfr <- scale(res.asfr, center=FALSE, scale=colSums(res.asfr))
	if(start.phase3 <= lyears) res.asfr <- update.by.mac(res.asfr, max(1, start.phase3-tobs))
	return(res.asfr)
}

.get.par.from.inputs <- function(par, inputs, country) {
	if(is.null(inputs[[par]])) return (NULL)
	idx <- inputs[[par]][,'country_code'] == country
	if(sum(idx)==0) return (NULL)
	res <- inputs[[par]][idx,,drop=FALSE]
	return (as.matrix(res[, !is.element(colnames(res), c('country_code', 'age')),drop=FALSE]))
}


get.country.inputs <- function(country, inputs, nr.traj, country.name) {
	inpc <- list()
	obs <- list()
	for(par in c('POPm0', 'POPf0', 'MXm', 'MXf', 'MXpattern', 'SRB',
				'PASFR', 'PASFRpattern', 'MIGtype', 'MIGm', 'MIGf', 'MXm.pred', 'MXf.pred')) {
		inpc[[par]] <- .get.par.from.inputs(par, inputs, country)
		obs[[par]] <- .get.par.from.inputs(par, inputs$observed, country)
	}
	inpc[['MIGBaseYear']] <- inpc[['MIGtype']][,'ProjFirstYear']
	inpc[['MIGtype']] <- inpc[['MIGtype']][,'MigCode']
	# generate sex and age-specific migration if needed
	if((!is.null(inpc[['MIGm']]) & all(is.na(inpc[['MIGm']]))) || (!is.null(inpc[['MIGf']]) & all(is.na(inpc[['MIGf']])))) {
		mig.recon <- age.specific.migration(wpp.year=inputs$wpp.year, countries=country, verbose=FALSE)
		mig.pair <- list(MIGm="male", MIGf="female")
		for(what.mig in names(mig.pair)) {
			if(!is.null(inpc[[what.mig]]) & all(is.na(inpc[[what.mig]]))) {
				# extact predicted migration
				cols <- intersect(colnames(mig.recon[[mig.pair[[what.mig]]]]), colnames(inpc[[what.mig]]))
				inpc[[what.mig]][,cols] <- as.matrix(mig.recon[[mig.pair[[what.mig]]]][,cols])
				rownames(inpc[[what.mig]]) <- rownames(mig.recon[[mig.pair[[what.mig]]]])
				# extact observed migration
				if(!is.null(obs[[what.mig]])) {
					cols <- intersect(colnames(mig.recon[[mig.pair[[what.mig]]]]), colnames(obs[[what.mig]]))
					obs[[what.mig]][,cols] <- as.matrix(mig.recon[[mig.pair[[what.mig]]]][,cols])
					rownames(obs[[what.mig]]) <- rownames(mig.recon[[mig.pair[[what.mig]]]])
				}
			}
		}
	}

	what.traj <- list(TFRpred='TFR', e0Mpred='male e0', e0Fpred='female e0')
	medians <- list()
	lyears <- length(inputs$proj.years)
	for(par in names(what.traj)) {
		if (!is.data.frame(inputs[[par]])) next
		cidx <- inputs[[par]][,'country_code'] == country 
		idx <- cidx & is.element(inputs[[par]][,'year'], inputs$proj.years)
		if(sum(idx) == 0) {
			warning('No ', what.traj[[par]], ' trajectories for ', country.name, 
				'. No population projection generated.')
			return(NULL)
		}
		df <- inputs[[par]][idx,-1]
		utrajs <- sort(unique(df$trajectory))
		ntrajs <- length(utrajs)
		if(ntrajs > 1){
			sorted.df <- data.frame(year=rep(inputs$proj.years, times=ntrajs), trajectory=rep(utrajs, each=lyears))
			# this is to get rows of the data frame in a particular order (align years and trajectories)
			df <- merge(sorted.df, df, sort=FALSE)
		}
		inpc[[par]] <- matrix(df[, 'value'], nrow=lyears)
		rownames(inpc[[par]]) <- as.character(inputs$proj.years)
		medians[[par]] <- apply(inpc[[par]], 1, quantile, 0.5, na.rm = TRUE)
		obsidx <- cidx & inputs[[par]][,'year'] < min(inputs$proj.years)
		if(sum(obsidx) > 0) {
			obs[[par]] <- inputs[[par]][obsidx,]
			obs.years <- unique(obs[[par]][, 'year'])
			obs[[par]] <- matrix(obs[[par]][, 'value'], nrow=length(obs.years))[,1]
			names(obs[[par]]) <- as.character(obs.years)
		}
	}
	e <- new.env()
	e$inputs <- inputs
	for(sex in c('M', 'F')) {
		par <- paste0('mig', sex, 'pred')
		if(is.null(inputs[[par]])) next
		inpc[[par]] <- .get.migration.traj(e, par, country)
		if(is.null(inpc[[par]])) next
		medians[[par]] <- apply(inpc[[par]], c(1,3), quantile, 0.5, na.rm = TRUE)
	}
	inpc$migMmedian <- medians$migMpred
	inpc$migFmedian <- medians$migFpred
	
	if(is.null(inpc$TFRpred)) {
		inpc$TFRpred <- get.tfr.trajectories(inputs$TFRpred, country)
		if(is.null(inpc$TFRpred)) {
			warning('No TFR trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)	
		}
		country.obj <- get.country.object(country, inputs$TFRpred$mcmc.set$meta)
		medians$TFRpred <- bayesTFR::get.median.from.prediction(inputs$TFRpred, country.obj$index, country.obj$code)[-1]
		obs.tfr <- bayesTFR:::get.tfr.reconstructed(inputs$TFRpred$tfr_matrix_reconstructed, inputs$TFRpred$mcmc.set$meta)
		obs$TFRpred <- obs.tfr[1:if(!is.null(inputs$TFRpred$present.year.index)) inputs$TFRpred$present.year.index else nrow(obs.tfr),country.obj$index]
	} 
	inpc$TFRhalfchild <- bayesTFR:::get.half.child.variant(median=medians$TFRpred, increment=c(0.25, 0.4, 0.5))
	if(!all(is.element(inputs$proj.years, colnames(inpc$TFRhalfchild))))
		stop('Mismatch in projection periods of TFR and the target projection years.')
	inpc$TFRhalfchild <- inpc$TFRhalfchild[,as.character(inputs$proj.years)]
	
	if(!inputs$fixed.mx){
	if(is.null(inpc$e0Mpred)) {
		inpc$e0Mpred <- get.e0.trajectories(inputs$e0Mpred, country)
		if(is.null(inpc$e0Mpred)) {
			warning('No male e0 trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)	
		}
		country.obj <- get.country.object(country, inputs$e0Mpred$mcmc.set$meta)
		medians$e0Mpred <- bayesTFR::get.median.from.prediction(inputs$e0Mpred, country.obj$index, country.obj$code)
		obs.e0M <- bayesLife:::get.e0.reconstructed(inputs$e0Mpred$e0.matrix.reconstructed, inputs$e0Mpred$mcmc.set$meta)
		obs$e0Mpred <- obs.e0M[1:if(!is.null(inputs$e0Mpred$present.year.index)) inputs$e0Mpred$present.year.index else nrow(obs.e0M),country.obj$index]
	}
	if(!all(is.element(inputs$proj.years, names(medians$e0Mpred))))
		stop('Mismatch in projection periods of male e0 and the target projection years.')
	inpc$e0Mmedian <- medians$e0Mpred[as.character(inputs$proj.years)]
	
	if(is.null(inpc$e0Fpred)) {
		inpc$e0Fpred <- get.e0.trajectories(inputs$e0Fpred, country)
		if(is.null(inpc$e0Fpred)) {
			warning('No female e0 trajectories for ', country.name, 
					'. No population projection generated.')
			return(NULL)
		}
		country.obj <- get.country.object(country, inputs$e0Fpred$mcmc.set$meta)
		medians$e0Fpred <- bayesTFR::get.median.from.prediction(inputs$e0Fpred, country.obj$index, country.obj$code)
		obs.e0F <- bayesLife:::get.e0.reconstructed(inputs$e0Fpred$e0.matrix.reconstructed, inputs$e0Fpred$mcmc.set$meta)
		obs$e0Fpred <- obs.e0F[1:if(!is.null(inputs$e0Fpred$present.year.index)) inputs$e0Fpred$present.year.index else nrow(obs.e0F),country.obj$index]

	}
	if(!all(is.element(inputs$proj.years, names(medians$e0Fpred))))
		stop('Mismatch in projection periods of female e0 and the target projection years.')
	inpc$e0Fmedian <- medians$e0Fpred[as.character(inputs$proj.years)]
	}
	# select trajectories
	indices <- list()
	nr.traj <- min(nr.traj, max(ncol(inpc$e0Mpred), ncol(inpc$e0Fpred), ncol(inpc$TFRpred), 
							    ncol(inpc$migMpred), ncol(inpc$migFpred)))
	for (par in c('TFRpred', 'e0Fpred', 'migMpred')) {
		if(is.null(inpc[[par]])) next
		traj.available <- ncol(inpc[[par]])
		if(traj.available == nr.traj) {
			indices[[par]] = 1:nr.traj
			next
		}
		if(traj.available > nr.traj) # select equidistantly
			indices[[par]] <- get.traj.index(nr.traj, inpc[[par]]) 
		else # re-sample
			indices[[par]] <- sample(1:traj.available, nr.traj, replace=TRUE)
	}
	# dependent trajectories
	dependencies <- list(e0Fpred='e0Mpred', migMpred='migFpred')
	for(par in names(dependencies)) {
		if(is.null(inpc[[dependencies[[par]]]])) {
			if(par == 'migMpred') {
				if((is.null(inpc[[par]]) || ncol(inpc[[par]])==1)) next
				stop('Female migration trajectories must match male migration.')
			}
			next
		}		
		if(par == 'migMpred' && (is.null(inpc[[par]]) || (!is.null(inpc[[par]]) && ncol(inpc[[par]])==1))) {
			indices[[dependencies[[par]]]] <- rep(1, nr.traj)
			next # using one trajectory, male is default
		}
		traj.available <- ncol(inpc[[par]])
		if (traj.available == ncol(inpc[[dependencies[[par]]]])) # same number of trajectories, therefore same indices
			indices[[dependencies[[par]]]] <- indices[[par]]
		else
			stop('Trajectories for female ', list(e0Fpred='life expectancy', migFpred='migration')[[par]], ' do not match the male ones: ',
				traj.available, ' <> ', ncol(inpc[[dependencies[[par]]]]), ' for ', country.name, '. No population projection generated.')
	}
	for (par in c('TFRpred', 'e0Mpred', 'e0Fpred')) { # these are matrices
		if(is.null(inpc[[par]])) next
		inpc[[par]] <- inpc[[par]][,indices[[par]], drop=FALSE]
		inpc[[par]] <- inpc[[par]][as.character(inputs$proj.years),, drop=FALSE]
	}
	inpc$mig.nr.traj <- 1
	for(par in c('migMpred', 'migFpred')) { # age-specific, thus 3-d arrays
		if(is.null(inpc[[par]])) next
		inpc[[par]] <- inpc[[par]][,indices[[par]], , drop=FALSE]
		inpc$mig.nr.traj <- length(indices[[par]])
	}
	inpc$observed <- obs
	inpc$trajectory.indices <- indices
	return(inpc)
}

get.traj.index <- function(nr.traj, traj, traj.dim=2, sample=FALSE) {
	traj.tot <- dim(traj)[traj.dim]
	if(nr.traj >= traj.tot) 
		return(if(sample) sample(1:traj.tot, nr.traj, replace=TRUE) else 1:traj.tot)
	return(if(sample) sample(1:traj.tot, nr.traj) else as.integer(seq(1, traj.tot, length=nr.traj)))
}

rotateLC <- function(e0, bx, bux, axM, axF, e0u=102, p=0.5) {
	# Rotation from Li, Lee, Gerland, 2013
	lmin <- -12
	lmax <- 0
	machine.max <- log(.Machine$double.xmax)
	machine.min <- log(.Machine$double.xmin)
	npred <- length(e0)
	kranges <- ax <- list()
	ku <- rep(NA, npred)
	for(sex in 1:2)
		kranges[[sex]] <- list(ku=ku, kl=ku)
	ax <- list(axM, axF)
	#if(!is.null(dim(bx)) && length(dim(bx))==2) { # aids country; bx is already a matrix
	#	Bxt <- bx
	#} else {
		wt <- (e0 - 80)/(e0u-80)
		wst <- (0.5*(1+(sin(pi/2*(2*wt-1)))))^p
		Bxt <- matrix(NA, nrow=length(bx), ncol=npred)
		for(t in 1:npred) {
			Bxt[,t] <- switch(cut(e0[t], c(0, 80, 102, 9999), labels=FALSE, right=FALSE),
					bx, 
					(1-wst[t])*bx + wst[t]*bux,
					bux)
		}
	#}
	bx.lim=rbind(apply(Bxt, 2, function(x) min(x[x>0])), apply(Bxt, 2, function(x) max(x[x>0])))
	for(sex in 1:2) {
		kranges[[sex]]$kl <- pmin((lmax - apply(ax[[sex]], 2, min))/bx.lim[2,], machine.max)
		kranges[[sex]]$ku <- pmax((lmin - apply(ax[[sex]], 2, max))/bx.lim[1,], machine.min)
	}
	return(list(bx=Bxt, kranges=kranges))
}

modifiedLC <- function (npred, mxKan, eopm, eopf, verbose=FALSE, debug=FALSE) {
	eop  <- list(eopm, eopf)
    sr <- LLm <- list(matrix(0, nrow=27, ncol=npred), matrix(0, nrow=27, ncol=npred))
    Mx <- lx <- list(matrix(0, nrow=28, ncol=npred), matrix(0, nrow=28, ncol=npred))
    #Get the projected kt from eo, and make projection of Mx
    nproj <- npred
    rotKan <- rotateLC(0.5*(eop[[1]]+eop[[2]]), mxKan$bx, mxKan$bux, mxKan$male$axt, mxKan$female$axt)
    for (mxYKan in list(mxKan$female, mxKan$male)) { # iterate over male and female
    	#print(c('sex: ', mxYKan$sex))
    	#stop('')	 
    	res <- .C("LC", as.integer(nproj), as.integer(mxYKan$sex), as.numeric(mxYKan$axt), 
    		as.numeric(rotKan$bx), as.numeric(eop[[mxYKan$sex]]), 
    		Kl=as.numeric(rotKan$kranges[[mxYKan$sex]]$kl), Ku=as.numeric(rotKan$kranges[[mxYKan$sex]]$ku), 
			constrain=as.integer(mxYKan$sex == 1), 
			#constrain=as.integer(0),
			FMx=as.numeric(Mx[[2]]), FEop=as.numeric(eop[[2]]),
			LLm = as.numeric(LLm[[mxYKan$sex]]), Sr=as.numeric(sr[[mxYKan$sex]]), 
			lx=as.numeric(lx[[mxYKan$sex]]), Mx=as.numeric(Mx[[mxYKan$sex]]))
		sr[[mxYKan$sex]] <- matrix(res$Sr, nrow=27)
		LLm[[mxYKan$sex]] <- matrix(res$LLm, nrow=27)
		Mx[[mxYKan$sex]] <- matrix(res$Mx, nrow=28)
		lx[[mxYKan$sex]] <- matrix(res$lx, nrow=28)
    }
	return(list(sr=sr, LLm=LLm, mx=Mx, lx=lx))    
}

survival.fromLT <- function (npred, mxKan, verbose=FALSE, debug=FALSE) {
    sr <- LLm <- list(matrix(0, nrow=27, ncol=npred), matrix(0, nrow=27, ncol=npred))
    Mx <- lx <- list(matrix(0, nrow=28, ncol=npred), matrix(0, nrow=28, ncol=npred))
    sx <- rep(0, 27)
    for (mxYKan in list(mxKan$female, mxKan$male)) { # iterate over male and female
    	Mx[[mxYKan$sex]] <- mxYKan$mx
    	sex <- c('Male', 'Female')[mxYKan$sex]
    	for(time in 1:npred) {
			res <- LifeTableMx(mxYKan$mx[,time], sex=sex)
			LLm[[mxYKan$sex]][,time] <- c(sum(res$Lx[1:2]), +res$Lx[3:nrow(res)]) # collapse first two age groups
			res.sr <- .C("get_sx27", as.numeric(LLm[[mxYKan$sex]][,time]), sx=sx)
			sr[[mxYKan$sex]][,time] <- res.sr$sx			
			lx[[mxYKan$sex]][,time] <- res$lx
		}
    }
    #stop('')
	return(list(sr=sr, LLm=LLm, mx=Mx, lx=lx))    
}

runKannisto <- function(inputs, start.year, ...) {
	# extend mx, get LC ax,bx,k1
	Kan <- KannistoAxBx.joint(inputs$MXm, inputs$MXf, start.year=start.year, mx.pattern=inputs$MXpattern, ...)
	mxMKan <- c(Kan$male, sex=1)
	mxFKan <- c(Kan$female, sex=2)
	bux <- NULL
	#if(is.null(mxMKan$bxt)) {
		bx <- 0.5 * (mxMKan$bx + mxFKan$bx)
		# ultimate bx (Li, Lee, Gerland 2013)
    	bux <- bx
    	avg15.65 <- mean(bux[5:14])
    	bux[1:14] <- avg15.65
    	bux[15:28] <- bux[15:28] * (bux[14]/bux[15]) # adjust so that b(70)=b(65)
    	bux <- bux/sum(bux) # must sum to 1
    #} else # aids country, bxt is a matrix
    #	bx <- 0.5 * (mxMKan$bxt + mxFKan$bxt)
	return(list(male=mxMKan, female=mxFKan, bx=bx, bux=bux))
}

runKannisto.noLC <- function(inputs) {
	# extend mx
	Kan <- KannistoAxBx.joint(inputs$MXm.pred, inputs$MXf.pred, compute.AxBx=FALSE)
	mxMKan <- c(Kan$male, sex=1)
	mxFKan <- c(Kan$female, sex=2)
	return(list(male=mxMKan, female=mxFKan))
}


KannistoAxBx.joint <- function(male.mx, female.mx, start.year=1950, mx.pattern=NULL, ax.from.latest.periods=99, npred=19, 
								joint=TRUE, compute.AxBx=TRUE)  {
	# Extending mx to age 130 using Kannisto model and mx 80-99, OLS
	finish.bx <- function(bx) {
			negbx <- which(bx <= 0)
			lnegbx <- length(negbx)
			if(lnegbx > 0 && negbx[1] == 1) {
				bx[1] <- 0
				negbx <- if(lnegbx > 1) negbx[2:lnegbx] else c()
				lnegbx <- length(negbx)
			}
			while(lnegbx > 0) {
				bx[negbx] <- 0.5 * bx[negbx-1]
				negbx <- which(bx < 0)
				lnegbx <- length(negbx)
			}      
			for (i in 1:27) { 
				if (bx[29 - i] == 0) bx[29 - i] <- bx[29 - i - 1]
			}
			bx <- bx/sum(bx) # must sum to 1
			return(bx)
		}

	Mxe.m <- as.matrix(male.mx)
	Mxe.m <- rbind(Mxe.m, 
			matrix(NA, nrow=28-nrow(male.mx), ncol=ncol(male.mx)))			
	Mxe.f <- as.matrix(female.mx)
	Mxe.f <- rbind(Mxe.f, 
			matrix(NA, nrow=28-nrow(female.mx), ncol=ncol(female.mx)))
	ne <- ncol(Mxe.m)
	k <- 22:28
	npoints <- 4
	age.group <- (21-npoints+1):21
	data <- data.frame(sex=c(rep(1,npoints), rep(0,npoints)), age=c(age.group, age.group))
	mxc <- rbind(male.mx[age.group, 1:ne, drop=FALSE], female.mx[age.group, 1:ne, drop=FALSE])
	logit.mxc <- log(mxc) - log(1-mxc)
	if(joint) {
		for(j in 1:ne) {		
			data$lmx <- logit.mxc[,j]
			if(all(is.na(data$lmx))) next
			fit <- lm(lmx ~ sex + age, data=data)
			coefs <- coefficients(fit)
			aam.female <- exp(coefs[1]) # intercept
			aam.male <- exp(coefs[1] + coefs['sex'])
			bbm <- coefs['age']
			# Ages 100-105, ..., 130+
			expterm.m <- aam.male * exp(bbm * k)	
			Mxe.m[k, j] =  expterm.m / (1 + expterm.m)
			expterm.f <- aam.female * exp(bbm * k)	
			Mxe.f[k, j] =  expterm.f / (1 + expterm.f)
		}
	} else {
		h <- 100 + 5 * (0:6)
		hminus80 <- h - 80
		lmxr <- list(M=log(male.mx[18:21, 1:ne] / (1 - male.mx[18:21, 1:ne])),
					F=log(female.mx[18:21, 1:ne] / (1 - female.mx[18:21, 1:ne])))
		res <- list(M=Mxe.m, F=Mxe.f)
		for(sex in c('M', 'F')) {
			Xm1 <- apply(lmxr[[sex]], 2, sum)
			Xm2 <- apply((5 * (1:4) - 5) * lmxr[[sex]], 2, sum)
			for(j in 1:ne) {		
				aam <- exp((350 * Xm1[j] - 30 * Xm2[j]) / 500)
				bbm <- (Xm1[j] - 4 * log(aam)) / 30
				expterm <- aam * exp(bbm * (hminus80))
				res[[sex]][k, j] =  expterm / (1 + expterm)
			}
		}
		Mxe.m <- res$M
		Mxe.f <- res$F
	}
	result <- list(male=list(mx=Mxe.m), female=list(mx=Mxe.f))
	if(!compute.AxBx) return(result)
	#Get Lee-Cater Ax and Bx
	years <- substr(colnames(male.mx),1,4)
	ns <- which(years == start.year)
	if(length(ns)==0) stop('start.year must be between ', years[1], ' and ', years[ne])
    model.bx <- !is.null(mx.pattern) && "AgeMortalityType" %in% colnames(mx.pattern) && mx.pattern[,"AgeMortalityType"] == "Model life tables"
    avg.ax <- !is.null(mx.pattern) && "LatestAgeMortalityPattern" %in% colnames(mx.pattern) && mx.pattern[,"LatestAgeMortalityPattern"] == 0
    smooth.ax <-  !is.null(mx.pattern) && !avg.ax && "SmoothLatestAgeMortalityPattern" %in% colnames(mx.pattern) && mx.pattern[,'SmoothLatestAgeMortalityPattern'] == 1
    is.aids.country <- !is.null(mx.pattern) && "WPPAIDS" %in% colnames(mx.pattern) && mx.pattern[,"WPPAIDS"] == 1
    if(is.aids.country) {
    	avg.ax <- FALSE
    	smooth.ax <- TRUE
    	aids.idx <- which(years < 1985)
    	aids.npred <- min((2100-(as.integer(years[ne])+5))/5, npred)
    }
    #avg.ax <- TRUE
    if(!avg.ax) ax.from.latest.periods <- 1
    mlt.bx <- NULL
    if(model.bx) {
    	bx.env <- new.env()
    	data(MLTbx, envir = bx.env)
    	bx.pattern <- if ("AgeMortalityPattern" %in% colnames(mx.pattern)) mx.pattern[,"AgeMortalityPattern"] else "UN General"
    	mlt.bx <- as.numeric(bx.env$MLTbx[bx.pattern,])
    }
    
    for(sex in c('male', 'female')) {
    	lMxe <- log(result[[sex]]$mx)
    	this.ns <- if(any(is.na(lMxe[,ns:ne]))) ns + sum(apply(lMxe[,ns:ne], 2, function(z) all(is.na(z))))
    				else ns
    	ax.ns <- max(ne-ax.from.latest.periods+1, this.ns)
    	x1 <- apply(lMxe[,ax.ns:ne, drop=FALSE], 1, sum, na.rm=TRUE)
    	ax <- x1 / (ne - ax.ns + 1)
    	if(smooth.ax) {
    		ax.sm <- smooth.spline(ax[1:21], df=11)$y
    		ax[2:21] <- ax.sm[2:21] # keep value the first age group
    	}
		kt <- rep(NA, ne)
		kt[this.ns:ne] = apply(lMxe[,this.ns:ne, drop=FALSE], 2, sum) - sum(ax)
		axt <- matrix(ax, nrow=28, ncol=npred)
		if(is.aids.country) {
			ax.end <- apply(lMxe[,aids.idx, drop=FALSE], 1, sum, na.rm=TRUE)/length(aids.idx)
			ax.end.sm <- smooth.spline(ax.end[1:21], df=11)$y
    		ax.end[2:21] <- ax.end.sm[2:21] # keep value the first age group
			for (i in 1:28) { # linear interpolation to the average ax ending in 2050; after that the avg ax is used
				axt[i,1:aids.npred] <- approx(c(1,aids.npred), c(ax[i], ax.end[i]), xout=1:aids.npred)$y
				if(aids.npred < npred)
					axt[i,(aids.npred+1):npred] <- ax.end[i]	
			}
		}
		bx <- mlt.bx
		#bxt <- NULL
    	if(!model.bx) {
			x2 <- sum(kt[this.ns:ne]*kt[this.ns:ne])
			x1 <- rep(NA, nrow(lMxe))
			for (i in 1:nrow(lMxe)) 
				x1[i] <- sum((lMxe[i,this.ns:ne]-ax[i])*kt[this.ns:ne])
			bx <- x1/x2
			bx <- finish.bx(bx)
			# if(is.aids.country) {
				# bxt <- matrix(bx, nrow=28, ncol=npred)
				# slMxe.aids <- apply(lMxe[,aids.idx, drop=FALSE], 2, sum)
				# x1.t <- rep(NA, nrow(lMxe))
				# for(t in 2:aids.npred) {
					# kt.t = slMxe.aids - sum(axt[,t])
					# x2.t <- sum(kt.t*kt.t)					
					# for (i in 1:nrow(lMxe)) 
						# x1.t[i] <- sum((lMxe[i,aids.idx]-axt[i,t])*kt.t)
					# bxt.t <- x1.t/x2.t
					# bxt[,t] <- finish.bx(bxt.t)
				# }
				# for(t in (aids.npred+1):npred) bxt[,t] <- bxt[,aids.npred]
			# }
		}
		result[[sex]]$ax <- ax
		result[[sex]]$axt <- axt
		result[[sex]]$bx <- bx
		#result[[sex]]$bxt <- bxt
		#result[[sex]]$k0 <- kt[ne]
		#result[[sex]]$d1 <- (kt[ne] - kt[this.ns]) / (ne - this.ns + 1)		
	}
	return(result)
}


StoPopProj <- function(npred, inputs, LT, asfr, mig.pred=NULL, mig.type=NULL, country.name=NULL, keep.vital.events=FALSE) {
	popm <- popf <- matrix(0, nrow=27, ncol=npred+1)
	popm[,1] <- c(inputs$POPm0, rep(0, 6))
	popf[,1] <- c(inputs$POPf0, rep(0, 6))
	totp <- c(sum(popm[,1]+popf[,1]), rep(0, npred))
	btageM <- btageF <- matrix(0, nrow=7, ncol=npred) # births by age of mother and sex of child
	deathsM <- deathsF <- matrix(0, nrow=27, ncol=npred)
	nproj <- npred
	migM <- as.matrix(if(!is.null(mig.pred[['M']])) mig.pred[['M']] else inputs[['MIGm']])
	migF <- as.matrix(if(!is.null(mig.pred[['F']])) mig.pred[['F']] else inputs[['MIGf']])
	#migM <- as.matrix(if(!is.null(mig.pred[['M']])) apply(mig.pred[['M']], 2, '/', popm[1:21,1]) else inputs[['MIGm']])
	#migF <- as.matrix(if(!is.null(mig.pred[['F']])) apply(mig.pred[['F']], 2, '/', popf[1:21,1]) else inputs[['MIGf']])
	#stop('')
	res <- .C("TotalPopProj", as.integer(nproj), as.numeric(migM), as.numeric(migF), nrow(migM), ncol(migM), as.integer(mig.type),
		srm=LT$sr[[1]], srf=LT$sr[[2]], asfr=as.numeric(as.matrix(asfr)), 
		srb=as.numeric(as.matrix(inputs$SRB)), 
		mxm=LT$mx[[1]], mxf=LT$mx[[2]],
		popm=popm, popf=popf, totp=totp,
		btagem=as.numeric(btageM), btagef=as.numeric(btageF), 
		deathsm=as.numeric(deathsM), deathsf=as.numeric(deathsF)
		)
	
 	if(any(res$popm < 0)) warnings('Negative male population for ', country.name)
	if(any(res$popf < 0)) warnings('Negative female population for ', country.name)

	vital.events <- list()
	if(keep.vital.events) {
		vital.events$mbt <- res$btagem
		vital.events$fbt <- res$btagef
		vital.events$mdeaths <- res$deathsm
		vital.events$fdeaths <- res$deathsf
	}
	return(c(list(totpop=res$totp, mpop=res$popm, fpop=res$popf), vital.events))
}

compute.observedVE <- function(inputs, pop.matrix, mig.type, mxKan, country.code, estim.years) {
	obs <- inputs$observed
	if(is.null(obs$PASFR)) return(NULL)
	npasfr <- nrow(obs$PASFR)
	nest <- min(length(obs$TFRpred), ncol(obs$PASFR), sum(!is.na(obs$MIGm[1,])), length(estim.years))
	estim.years <- estim.years[(length(estim.years)-nest+1):length(estim.years)]
	pasfr <- obs$PASFR[,(ncol(obs$PASFR)-nest+1):ncol(obs$PASFR), drop=FALSE]
	tfr <- obs$TFRpred[(length(obs$TFRpred)-nest+1):length(obs$TFRpred)]
	mig.data <- list(as.matrix(obs$MIGm[,(ncol(obs$MIGm)-nest+1):ncol(obs$MIGm)]), 
					as.matrix(obs$MIGf[,(ncol(obs$MIGf)-nest+1):ncol(obs$MIGf)]))
	asfr <- pasfr/100.
	for(i in 1:npasfr) asfr[i,] <- tfr * asfr[i,]
	
	pop <- D10 <- list()
	nmx <- ncol(inputs$MXm)
	mx <-  list(inputs$MXm[,(nmx-nest+1):nmx, drop=FALSE], inputs$MXf[,(nmx-nest+1):nmx, drop=FALSE])
	srb <- obs$SRB[(length(obs$SRB)-nest+1):length(obs$SRB)]
	srb.ratio <- srb / (1 + srb)
	sr <- list(matrix(0, nrow=21, ncol=nest), matrix(0, nrow=21, ncol=nest))
	births <- deaths <- list()
	for(sex in 1:2) {
		tpop <- get.pop.observed.with.age(NULL, country.code, sex=c('male', 'female')[sex], 
						data=pop.matrix)
		pop[[sex]] <- tpop$data[tpop$age.idx,(ncol(tpop$data)-nest):ncol(tpop$data)]
	}
	bt <- (pop[[2]][4:10,-1] + pop[[2]][4:10,-ncol(pop[[2]])]) * asfr * 0.5
	births[[1]] <- bt * srb.ratio
	births[[2]] <- bt - births[[1]]
	for(sex in 1:2) {		
		#sr[[sex]] <- get.survival(abind(mx[[sex]],along=3), sex=c("M","F")[sex], age05=c(TRUE, TRUE, FALSE))[,,1]
		sr[[sex]] <- get.survival(abind(mx[[sex]],along=3), sex=c("M","F")[sex])[1:21,,1]
		deaths[[sex]] <- matrix(0, nrow=21, ncol=nest)
		res <- .C("get_deaths_from_sr", as.numeric(sr[[sex]]), nest, as.numeric(as.matrix(pop[[sex]])), 
					as.numeric(mig.data[[sex]]), mig.type,
					as.numeric(colSums(as.matrix(births[[sex]]))), Deaths=as.numeric(deaths[[sex]]))
		#stop('')
		deaths[[sex]] <- matrix(res$Deaths, nrow=21)
		colnames(deaths[[sex]]) <- estim.years
		rownames(deaths[[sex]]) <- rownames(pop[[sex]])
		colnames(births[[sex]]) <- estim.years
	}	
	colnames(asfr) <- estim.years
	rownames(asfr) <- rownames(births[[1]])
	res <- list(btm=births[[1]], btf=births[[2]], 
				deathsm=deaths[[1]], deathsf=deaths[[2]], asfert=asfr, pasfert=pasfr, mxm=mx[[1]], mxf=mx[[2]])
	for(par in names(res))
		res[[par]] <- abind(res[[par]], NULL, along=3) # add trajectory dimension
	return(res)
}

write.pop.projection.summary <- function(pop.pred, what=NULL, expression=NULL, output.dir=NULL, ...) {
	pred <- as.environment(pop.pred) # to avoid lots of copying and allowing attaching cache
	if (is.null(output.dir)) output.dir <- pop.output.directory(pred)
	if(!file.exists(output.dir)) dir.create(output.dir, recursive=TRUE)
	all.what <- c('pop', 'popsex', 'popsexage', 'popage', 'births', 'birthssex', 'birthsage', 'birthssexage', 
			'deaths', 'deathssex', 'deathsage', 'deathssexage', 'srsexage', 'fertility', 'fertilityage', 'pfertilityage')
	what <- if(is.null(what)) all.what else match.arg(what, all.what, several.ok=TRUE)
	params <- list()
	if(!is.null(expression)) {
		what <- 'expression'
		params <- list(expression=expression)
	}
	for(summary.type in what) 
		do.call(paste0('write.', summary.type), c(list(pred, output.dir=output.dir), params, ...))
}

write.pop <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, ...)
	
write.popsex <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, what.log='population', ...)
	
write.popsexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, what.log='population', ...)
	
write.popage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, what.log='population', ...)
	
write.births <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='births', 
			file.suffix='births', what.log='total births', ...)
	
write.birthssex <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, vital.event='births', 
			file.suffix='births', what.log='births', ...)

write.birthsage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='births', 
			file.suffix='births', what.log='births', ...)

write.birthssexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='births', 
			file.suffix='births', what.log='births', ...)

write.deaths <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='deaths', 
			file.suffix='deaths', what.log='total deaths', ...)
	
write.deathssex <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=FALSE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths', ...)

write.deathsage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths', ...)
	
write.deathssexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='deaths', 
			file.suffix='deaths', what.log='deaths', ...)
			
write.srsexage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=TRUE, byage=TRUE, vital.event='survival', 
			file.suffix='sr', what.log='survival ratio', digits=litem('digits', list(...), 4))
			
write.fertility <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=FALSE, vital.event='fertility', 
			file.suffix='asfr', what.log='fertility rate', digits=litem('digits', list(...), 4))
			
write.fertilityage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='fertility', 
			file.suffix='asfr', what.log='fertility rate', digits=litem('digits', list(...), 4))
			
write.pfertilityage <- function(pop.pred, output.dir, ...) 
	.write.pop(pop.pred, output.dir=output.dir, bysex=FALSE, byage=TRUE, vital.event='pasfr', 
			file.suffix='pasfr', what.log='percent fertility rate', digits=litem('digits', list(...), 4))
	
write.expression <- function(pop.pred, expression, output.dir, file.suffix='expression', 
							expression.label=expression, include.observed=FALSE, digits=NULL, 
							adjust=FALSE, adj.to.file=NULL, end.time.only=FALSE) {
	cat('Creating summary file for expression ', expression, ' ...\n')
	header <- list(country.name='country_name',  country.code='country_code', variant='variant')
	variant.names <- c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95')
	nr.var <- length(variant.names)
	pred.period <- get.pop.prediction.periods(pop.pred, end.time.only=end.time.only)
	nr.proj <- length(pred.period)
	nr.obs <- 0
	if(include.observed) {
		obs.period <- get.pop.observed.periods(pop.pred, end.time.only=end.time.only)
		nr.obs <- length(obs.period)-1
		obs.period <- obs.period[-(nr.obs+1)] # remove the last one because the same as the first projection period
		for(iyear in 1:nr.obs) header[[paste0('year', iyear)]] <- obs.period[iyear]
	}
	for(iyear in 1:nr.proj) header[[paste0('year', iyear+nr.obs)]] <- pred.period[iyear]
	col.names <- grep('year', names(header), value=TRUE)
	result <- NULL
	if(include.observed) {
		res <- get.pop.observed.from.expression.all.countries(expression, pop.pred, time.index=1:nr.obs)
		#copy the same data into the variant rows 
		result <- matrix(NA, nrow=nrow(res)*5, ncol=ncol(res))
		for(i in 1:5) result[seq(i,by=5, length=nrow(res)),] <- res
	}
	if(adjust && is.null(pop.pred$adjust.env)) pop.pred$adjust.env <- new.env()
	for(iyear in 1:nr.proj) {	
		result <- cbind(result, as.vector(t(get.pop.from.expression.all.countries(expression, pop.pred, 
						quantiles=c(0.5, 0.1, 0.9, 0.025, 0.975), projection.index=iyear, adjust=adjust, adj.to.file=adj.to.file))))
	}
	if(!is.null(digits)) result <- round(result, digits)
	colnames(result) <- col.names
	result <- cbind(country.name=rep(as.character(pop.pred$countries$name), each=nr.var), 
					country.code=rep(pop.pred$countries$code, each=nr.var),
					variant=rep(variant.names, length(pop.pred$countries$code)), result)
	colnames(result)[colnames(result)==names(header)] <- header
	file <- file.path(output.dir, paste('projection_summary_', file.suffix, '.csv', sep=''))
	write(paste('#', expression.label), file)
	warn <- getOption('warn')
	options(warn=-1) # disable warning messages (it doesn't like that col.names is TRUE and append=TRUE)
	write.table(result, file=file, sep=',', row.names=FALSE, col.names=TRUE, append=TRUE,
				quote=which(is.element(colnames(result), c('country_name', 'variant'))))
	options(warn=warn)
	cat('Stored into: ', file, '\n')
}	
	
.write.pop <- function(pop.pred, output.dir, bysex=FALSE, byage=FALSE, vital.event=NULL, file.suffix='tpop', 
							what.log='total population', include.observed=FALSE, digits=0, adjust=FALSE, 
							end.time.only=FALSE) {
	cat('Creating summary file of ', what.log, ' ')
	if(bysex) cat('by sex ')
	if(byage) cat('by age ')
	cat ('...\n')
	header <- list(country.name='country_name',  country.code='country_code', variant='variant')
	if(bysex) header[['sex']] <- 'sex'
	if(byage) header[['age']] <- 'age'
	variant.names <- c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95')
	nr.var <- length(variant.names)
	if(missing(end.time.only)) end.time.only <- is.null(vital.event)
	pred.period <- get.pop.prediction.periods(pop.pred, end.time.only=end.time.only)
	#if(!is.null(vital.event)) pred.period <- pred.period[2:length(pred.period)]
	nr.obs <- 0
	if(include.observed) {
		obs.period <- get.pop.observed.periods(pop.pred, end.time.only=end.time.only)
		nr.obs <- length(obs.period)-1
		obs.period <- obs.period[-(nr.obs+1)] # remove the last one because the same as the first projection period
		for(iyear in 1:nr.obs) header[[paste0('year', iyear)]] <- obs.period[iyear]
	}
	nr.proj <- length(pred.period)
	for (i in 1:nr.proj) 
		header[[paste0('year', i+nr.obs)]] <- pred.period[i]
	col.names <- grep('year', names(header), value=TRUE)
	result <- NULL
	sex.index <- c(TRUE, FALSE, FALSE)
	if(bysex) sex.index <- !sex.index
	age.index <- c(TRUE, rep(FALSE, length(pop.pred$ages)))
	if(byage) age.index <- !age.index
	ages <- 1:length(pop.pred$ages)
	if(adjust && is.null(pop.pred$adjust.env)) pop.pred$adjust.env <- new.env()
	observed.data <- NULL
	for (country in 1:nrow(pop.pred$countries)) {
		country.obj <- get.country.object(country, country.table=pop.pred$countries, index=TRUE)
		for(sex in c('both', 'male', 'female')[sex.index]) {
			if(!is.null(vital.event)) {
			 	sum.over.ages <- age.index[1]
			 	if(include.observed) 
			 		observed <- get.popVE.trajectories.and.quantiles(pop.pred, country.obj$code, event=vital.event, 
										sex=sex, age='all', sum.over.ages=sum.over.ages, is.observed=TRUE)
				traj.and.quantiles <- get.popVE.trajectories.and.quantiles(pop.pred, country.obj$code, event=vital.event, 
										sex=sex, age='all', sum.over.ages=sum.over.ages)
				if(is.null(traj.and.quantiles$trajectories)) {
					warning('Problem with loading ', vital.event, '. Possibly no vital events stored during prediction.')
					return()
				}
				if(!sum.over.ages) { # This is because births have only subset of ages
					ages <- traj.and.quantiles$age.idx.raw
					age.index <- age.index[1:(length(ages)+1)]
					subtract.from.age <- traj.and.quantiles$age.idx.raw[1]-traj.and.quantiles$age.idx[1]
				}
				
			}
			for(age in c('all', ages)[age.index]) {
				this.result <- cbind(
							country.name=rep(country.obj$name, nr.var), 
							country.code=rep(country.obj$code, nr.var),
							variant=variant.names)
				if(sex != 'both')
					this.result <- cbind(this.result, sex=rep(sex, nr.var))
				if(age != 'all') {
					age <- as.integer(age)
					this.result <- cbind(this.result, age=rep(get.age.labels(pop.pred$ages)[age], nr.var))
				}
				if(is.null(vital.event)) {
					if(include.observed) 
						observed.data <- get.pop.observed(pop.pred, country.obj$code, sex=sex, age=age)
					quant <- get.pop.trajectories(pop.pred, country.obj$code, nr.traj=0, sex=sex, age=age, adjust=adjust)$quantiles
					traj <- NULL
					reload <- TRUE
				} else { # vital event
					quant <- traj.and.quantiles$quantiles
					traj <- traj.and.quantiles$trajectories
					if(include.observed)
						observed.data <- observed$trajectories[,,1]
					if(!sum.over.ages) {
						quant <- quant[age-subtract.from.age,,]
						traj <- traj[age-subtract.from.age,,]
						if(include.observed) {
							if (age-subtract.from.age > nrow(observed.data)) # because observed goes only up to 100+
								observed.data <- rep(0, ncol(observed.data))
							else
								observed.data <- observed.data[age-subtract.from.age,]
						}
					}
					reload <- FALSE
					#stop('')
				}
				proj.result <- round(rbind(
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, q=0.5, 
											trajectories=traj, reload=reload), 
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, pi=80, 
											trajectories=traj, reload=reload),
					get.pop.traj.quantiles(quant, pop.pred, country.obj$index, country.obj$code, pi=95, 
											trajectories=traj, reload=reload)),
					digits)
				if(!is.null(observed.data)) {
					# put it into the same shape as proj.result minus the last observed
					observed.data <- round(rbind(observed.data, NULL), digits)
					observed.data <- observed.data[rep(1, nrow(proj.result)), -ncol(observed.data)]
					proj.result <- cbind(observed.data, proj.result)
				}
				colnames(proj.result) <- col.names
				#stop('')
				this.result <- cbind(this.result, proj.result)
				result <- rbind(result, this.result)
			}
		}
	}
	colnames(result)[colnames(result)==names(header)] <- header
	suffix <- paste(file.suffix, paste(c('sex', 'age')[c(bysex, byage)], collapse=''), sep='')
	file <- paste('projection_summary_', suffix, '.csv', sep='')
	write.table(result, file=file.path(output.dir, file), sep=',', row.names=FALSE, col.names=TRUE, 
				quote=which(is.element(colnames(result), c('country_name', 'variant', 'sex', 'age'))))
	cat('Stored into: ', file.path(output.dir, file), '\n')
}

LifeTableMxCol <- function(mx, colname=c('Lx', 'lx', 'qx', 'mx', 'dx', 'Tx', 'sx', 'ex', 'ax'), ...){
	colname <- match.arg(colname)
	if(is.null(dim(mx))) return(.doLifeTableMxCol(mx, colname, ...))
	return(apply(mx, 2, .doLifeTableMxCol, colname=colname, ...))
}

.collapse.sx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	# sx does not need to be collapsed, as all elements refer to a five-year age group
	# for collapsed format remove last age group so to assure same length as all other columns after collapsing	
	#sx.start <- c(LT$sx[1:2], (LT$Lx[1] + LT$Lx[2])/5)[age05]
	#return(c(sx.start, LT$sx[-(1:2)]))
	return(switch(sum(age05)+1,
				LT$sx[c(-1, -length(LT$sx))], # unuseful case of all values FALSE
				LT$sx[-length(LT$sx)], # one value TRUE (collapsed life table)
				LT$sx, # two values TRUE (abridged life table)
				c(LT$sx[1], LT$sx) # unuseful case of all values TRUE
				))
}

.collapse.dx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	dx.start <- c(LT$dx[1:2], LT$dx[1] + LT$dx[2])[age05]
	return(c(dx.start, LT$dx[-(1:2)]))
}

.collapse.ax <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	ax.start <- c(LT$ax[1:2], ((LT$Lx[1]+LT$Lx[2]-5*LT$lx[3])/(LT$dx[1]+LT$dx[2])))[age05]
	return(c(ax.start, LT$ax[-(1:2)]))
}

.collapse.Tx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	Tx.start <- c(LT$Tx[c(1,2,1)])[age05]
	return(c(Tx.start, LT$Tx[-(1:2)]))
}

.collapse.ex <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	ex.start <- c(LT$ex[c(1,2,1)])[age05]
	return(c(ex.start, LT$ex[-(1:2)]))
}

.collapse.Lx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	Lx.start <- c(LT$Lx[1:2], LT$Lx[1] + LT$Lx[2])[age05]
	return(c(Lx.start, LT$Lx[-(1:2)]))
}

.collapse.lx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	lx.start <- LT$lx[c(1,2,1)][age05]
	return(c(lx.start, LT$lx[-(1:2)]))
}

.collapse.qx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	qx.start <- c(LT$qx[1:2], 1-(LT$lx[3]/LT$lx[1]))[age05]
	return(c(qx.start, LT$qx[-(1:2)]))
}

.collapse.mx <- function(LT, age05=c(FALSE, FALSE, TRUE)) {
	mx.start <- c(LT$mx[1:2], (LT$lx[1] - LT$lx[3])/(LT$Lx[1] + LT$Lx[2]))[age05]
	return(c(mx.start, LT$mx[-(1:2)]))
}

.doLifeTableMxCol <- function(mx, colname, age05=c(FALSE, FALSE, TRUE), ...) {
	# age05 determines the inclusion of ages 0-1, 1-4, 0-4
	LT <- LifeTableMx(mx, ...)
	if(all(age05==c(TRUE, TRUE, FALSE))) { # no collapsing
		result <- LT[,colname]
		names(result) <- rownames(LT)
	} else {
		result <- do.call(paste('.collapse', colname, sep='.'), list(LT, age05=age05))
		names(result) <- c(c('0-1', '1-4', '0-4')[age05], rownames(LT)[-(1:2)])
	}
	return(result)
}


LifeTableMx <- function(mx, sex=c('Male', 'Female'), include01=TRUE){
	# The first two elements of mx must correspond to 0-1 and 1-4. 
	# If include01 is FALSE, the first two age groups of the results are collapsed to 0-5
	sex <- match.arg(sex)
	sex <- list(Male=1, Female=2)[[sex]]
	nage <- length(mx)
	Lx <- lx <- qx <- Tx <- sx <- dx <- ax <- rep(0, nage)
	nagem1 <- nage-1
	resage <- c(if(include01) c(0,1) else 0, seq(5, by=5, length=nage-2))
	nresage <- length(resage)
	nas <- rep(NA,nresage)
	if(any(is.na(mx))) # there are NAs in mx
		return(data.frame(age=resage, mx=mx[1:nresage], qx=nas, lx=nas, dx=nas, Lx=nas, sx=nas, Tx=nas, ex=nas, ax=nas))
	
	LTC <- .C("LifeTable", as.integer(sex), as.integer(nagem1), as.numeric(mx), 
					Lx=Lx, lx=lx, qx=qx, ax=ax, Tx=Tx, sx=sx, dx=dx)
	LT <- data.frame(age=c(0,1, seq(5, by=5, length=nage-2)), 
					mx=mx, qx=LTC$qx, lx=LTC$lx, dx=LTC$dx, Lx=LTC$Lx,  sx=LTC$sx, Tx=LTC$Tx, ex=LTC$Tx/LTC$lx, ax=LTC$ax)
	LT$ax[nage] <- LT$ex[nage]
	if(!include01) { # collapse 0-1 and 1-5 into 0-5
		age05 <- c(FALSE, FALSE, TRUE)
		LTres <- data.frame(age=resage)
		for(colname in colnames(LT)[-1])
			LTres[[colname]] <- do.call(paste('.collapse', colname, sep='.'), list(LT, age05=age05))
		LT <- LTres
	}
	return(LT)
}

unblock.gtk.if.needed <- function(status, gui.opts=list()) {
	# This is to unblock the GUI, if the run is invoked from bayesDem
	# and pass info about its status
	# In such a case the gtk libraries are already loaded
	if(getOption('bDem.Poppred', default=FALSE)) {
		gui.opts$bDem.Poppred.status <- status
		unblock.gtk('bDem.Poppred', gui.opts)
	}
}

unblock.gtk <- function(...) bayesTFR:::unblock.gtk(...)

create.pop.cluster <- function(nr.nodes, ...) {
	cl <- makeCluster(nr.nodes, ...)
	lib.paths <- .libPaths()
	clusterExport(cl, c("lib.paths"), envir=environment())
	clusterEvalQ(cl, {.libPaths(lib.paths); library(bayesPop)})
	return(cl)
}

age.specific.migration <- function(wpp.year=2015, years=seq(1955, 2100, by=5), countries=NULL, smooth=TRUE, 
									write.to.disk=FALSE, directory=getwd(), file.prefix="migration", verbose=TRUE) {
	# Reconstruct sex- and age-specific net migration using a residual method using wpp data on population
	# and other available indicators. It is scaled to the total net migration for each country. 
	# It is not balanced over the world. Due to rounding issues, often it results in zig-zags over the ages,
	# therefore it is smoothed (in a double pass through the smoother).
	if(verbose) {
		status.text <- paste('Reconstructing sex- and age-specific migration from', paste0('wpp', wpp.year, ' '))
		cat('\n', status.text)
	}
	popm0 <- load.wpp.dataset("popM", wpp.year)
	popm0.num.cols <- grep('^[0-9]{4}$', colnames(popm0), value=TRUE) # values of year-columns
	popf0 <- load.wpp.dataset("popF", wpp.year)
	popf0.num.cols <- grep('^[0-9]{4}$', colnames(popf0), value=TRUE)
	popmproj <- load.wpp.dataset("popMprojMed", wpp.year)
	popmproj.num.cols <- grep('^[0-9]{4}$', colnames(popmproj), value=TRUE)
	popfproj <- load.wpp.dataset("popFprojMed", wpp.year)
	popfproj.num.cols <- grep('^[0-9]{4}$', colnames(popfproj), value=TRUE)
	sexrat <- load.wpp.dataset("sexRatio", wpp.year)
	sexrat.num.cols <- grep('^[0-9]{4}', colnames(sexrat), value=TRUE)
	mxm <- load.wpp.dataset("mxM", wpp.year)
	mxm.num.cols <- grep('^[0-9]{4}', colnames(mxm), value=TRUE)
	mxf <- load.wpp.dataset("mxF", wpp.year)
	mxf.num.cols <- grep('^[0-9]{4}', colnames(mxf), value=TRUE)
	mig <- load.wpp.dataset("migration", wpp.year)
	mig.num.cols <- grep('^[0-9]{4}', colnames(mig), value=TRUE)
	tfrproj <- .load.wpp.traj('tfr', wpp.year, median.only=TRUE)
	pasfr <- load.wpp.dataset("percentASFR", wpp.year)
	pasfr.num.cols <- grep('^[0-9]{4}', colnames(pasfr), value=TRUE)
	vwBase <- read.bayesPop.file(paste('vwBaseYear', wpp.year, '.txt', sep=''))[,c('country_code', 'MigCode')]
	pop.first.country <- popm0[popm0$country_code == mig$country_code[1],]
	max.ages <- nrow(pop.first.country)
	ages <- 1:max.ages
	age.labels <- get.age.labels(ages, age.is.index=TRUE, last.open=TRUE)
	years.periods <- paste(years-5, years, sep="-")
	lyears <- length(years)
	if(is.null(countries)) {
		countries <- mig$country_code
		# filter out non-countries
		locs <- bayesTFR:::load.bdem.dataset('UNlocations', wpp.year, envir=globalenv())
		countries <- countries[countries %in% locs[locs$location_type==4, "country_code"]]
	} else mig <- mig[which(mig$country_code %in% countries),]
	all.migM <- all.migF <- NULL
	lcountries <- length(countries)
	for(icountry in 1:lcountries) {
		if(verbose && interactive()) cat('\r', status.text, round(icountry/lcountries*100), '%')
		country <- countries[icountry]
		country.name <- as.character(mig$name[icountry])
		# filter country data
		popm.obs <- popm0[popm0$country_code==country, popm0.num.cols]
		popf.obs <- popf0[popf0$country_code==country, popf0.num.cols]
		pop1m <- cbind(popm.obs, popmproj[popmproj$country_code==country, popmproj.num.cols])
		pop1f <- cbind(popf.obs, popfproj[popfproj$country_code==country, popfproj.num.cols])
		tfra <- tfrproj[tfrproj$country_code==country,]
		asfr <- pasfr[pasfr$country_code==country, pasfr.num.cols]
		sr <- sexrat[sexrat$country_code==country, sexrat.num.cols] 
		mortM <- mxm[mxm$country_code==country, mxm.num.cols]
		mortF <- mxf[mxf$country_code==country, mxf.num.cols]
		totmig <- mig[mig$country_code==country, mig.num.cols]
		mtype <- vwBase[vwBase$country_code==country,'MigCode']
		if(length(mtype)==0) mtype <- 9
		this.all.migM <- this.all.migF <- data.frame(
				country_code=rep(country, max.ages), name=rep(country.name, max.ages), age=age.labels)
		for(iyear in 1:lyears) {
			year <- years[iyear]
			year.col <- years.periods[iyear]
			year.char <- as.character(year)
			pop0m <- pop1m[,as.character(year-5)]
			pop0f <- pop1f[,as.character(year-5)]
			mortMy <- mortM[,year.col]
			mortFy <- mortF[,year.col]
			sxm <- get.survival(matrix(mortMy, ncol=1), sex="Male")[,1,1]
       		sxf <- get.survival(matrix(mortFy, ncol=1), sex="Female")[,1,1]
			totmigy <- round(totmig[,year.col],3)
			if(totmigy == 0) netmigM <- netmigF <- rep(0, max.ages)
			else {			
				B2 <- sum((pop1f[4:10,year.char] + pop0f[4:10])/2 * tfra[tfra$year==year-2,'value'] * asfr[,year.col]/100)
				netmigM <- c(NA, pop1m[2:max.ages,year.char] - (pop0m[1:(max.ages-1)] * sxm[2:max.ages]))
				netmigF <- c(NA, pop1f[2:max.ages,year.char] - (pop0f[1:(max.ages-1)] * sxf[2:max.ages]))
				B2m <- B2 * sr[,year.col]/(1+sr[,year.col])
				netmigM[1] <- pop1m[1,year.char] - B2m * sxm[1]
				netmigF[1] <- pop1f[1,year.char] - (B2 - B2m) * sxf[1]
				migdata <- list(M=netmigM, F=netmigF)
				sxdata <- list(M=sxm, F=sxf)
				for(sex in c('M', 'F')) {
					if(mtype == 0) {					
					 	migdata[[sex]] <- 2*migdata[[sex]]
					 	for(i in 2:max.ages) {
					       migdata[[sex]][i] <- migdata[[sex]][i] - migdata[[sex]][i-1]*sxdata[[sex]][i]
					  	}
					}
					migdata[[sex]][18:21] <- 0
					if(smooth) { #smoothing					
						for(izig in 1:2) { # two passes of smoothing
							# are there significant zig-zags?						
							tops <- ((migdata[[sex]][2:(max.ages-1)] > migdata[[sex]][1:(max.ages-2)] & 
											migdata[[sex]][2:(max.ages-1)] > migdata[[sex]][3:max.ages]) | 
								(migdata[[sex]][2:(max.ages-1)] < migdata[[sex]][1:(max.ages-2)] & migdata[[sex]][2:(max.ages-1)] < migdata[[sex]][3:max.ages])) & 
								abs(diff(migdata[[sex]])[1:(max.ages-2)]/(totmigy/100)) > 0.05
							cs <- cumsum(c(TRUE, tops)) # consider first point as top
							if(any(cs[3:(max.ages-1)] > 2 & cs[3:(max.ages-1)] - cs[1:(max.ages-3)] > 1))  { # at least 3 neighboring tops
								migdata[[sex]] <- smooth.spline(migdata[[sex]], df=10)$y # smooth
								migdata[[sex]][18:21] <- 0
							} else break
						}
					}
				}
				# rescale
				netmigM <- migdata[['M']]
				netmigF <- migdata[['F']]
				s <- sum(netmigM + netmigF)
				netmigM <- netmigM/s * totmigy
				netmigF <- netmigF/s * totmigy
			}
			this.all.migM <- cbind(this.all.migM, netmigM)
			this.all.migF <- cbind(this.all.migF, netmigF)
		}
		colnames(this.all.migM) <- colnames(this.all.migF) <- c('country_code', 'name', 'age', years.periods)
		all.migM <- rbind(all.migM, this.all.migM)
		all.migF <- rbind(all.migF, this.all.migF)
	}
	if(write.to.disk) {
		write.table(all.migM, file=file.path(directory, paste0(file.prefix, "M.txt")), sep='\t', row.names=FALSE)
		write.table(all.migF, file=file.path(directory, paste0(file.prefix, "F.txt")), sep='\t', row.names=FALSE)
		if(verbose) cat('\nMigration files written into ', file.path(directory, paste0(file.prefix, "X.txt")))
	}
	if(verbose) cat('\n')
	return(invisible(list(male=all.migM, female=all.migF)))
}
