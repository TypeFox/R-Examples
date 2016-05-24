if(getRversion() >= "2.15.1") utils::globalVariables("loess_sd")
data(loess_sd, envir=environment())

e0.proj.le.SDPropToLoess<-function(x,l.start,kap,n.proj=11, p1=9, p2=9, const.var=FALSE){
  proj<-NULL
  proj[1]<-l.start
  for(a in 2:(n.proj+1)){
  	proj[a]<-proj[a-1]+g.dl6(x,proj[a-1], p1=p1, p2=p2)+rnorm(1,mean=0,
  				sd=(kap*if(const.var) 1 else loess.lookup(proj[a-1])))
  }
  return(proj[2:length(proj)])
}

e0.predict <- function(mcmc.set=NULL, end.year=2100, sim.dir=file.path(getwd(), 'bayesLife.output'),
                       replace.output=FALSE, predict.jmale = TRUE, nr.traj = NULL, thin=NULL, burnin=10000, 
                       use.diagnostics=FALSE, save.as.ascii=1000, start.year=NULL,
                       output.dir = NULL, low.memory=TRUE, seed=NULL, verbose=TRUE, ...){
	if(!is.null(mcmc.set)) {
		if (class(mcmc.set) != 'bayesLife.mcmc.set') {
			stop('Wrong type of mcmc.set. Must be of type bayesLife.mcmc.set.')
		}
	} else {                
		mcmc.set <- get.e0.mcmc(sim.dir, low.memory=low.memory, verbose=verbose)
	}
	if(!is.null(seed)) set.seed(seed)
		# Get argument settings from existing convergence diagnostics
	if(use.diagnostics) {
		diag.list <- get.e0.convergence.all(mcmc.set$meta$output.dir)
		ldiag <- length(diag.list)
		if (ldiag == 0) stop('There is no diagnostics available. Use manual settings of "nr.traj" or "thin".')
		use.nr.traj <- use.burnin <- rep(NA, ldiag)
		for(idiag in 1:ldiag) {
			if (bayesTFR::has.mcmc.converged(diag.list[[idiag]])) {
				use.nr.traj[idiag] <- diag.list[[idiag]]$use.nr.traj
				use.burnin[idiag] <- diag.list[[idiag]]$burnin
			}
		}
		if(all(is.na(use.nr.traj)))
			stop('There is no diagnostics indicating convergence of the MCMCs. Use manual settings of "nr.traj" or "thin".')
		# Try to select those that suggest nr.traj >= 2000 (take the minimum of those)
		traj.is.notna <- !is.na(use.nr.traj)
		larger2T <- traj.is.notna & use.nr.traj>=2000
		nr.traj.idx <- if(sum(larger2T)>0) (1:ldiag)[larger2T][which.min(use.nr.traj[larger2T])] 
						else (1:ldiag)[traj.is.notna][which.max(use.nr.traj[traj.is.notna])]
		nr.traj <- use.nr.traj[nr.traj.idx]
		burnin <- use.burnin[nr.traj.idx]
		if(verbose)
			cat('\nUsing convergence settings: nr.traj=', nr.traj, ', burnin=', burnin, '\n')
	}

	pred <- make.e0.prediction(mcmc.set, end.year=end.year,  
					replace.output=replace.output,  
					nr.traj=nr.traj, thin=thin, burnin=burnin, save.as.ascii=save.as.ascii, start.year=start.year,
					output.dir=output.dir, verbose=verbose)
	if(predict.jmale && mcmc.set$meta$sex == 'F')
		pred <- e0.jmale.predict(pred, ..., save.as.ascii=save.as.ascii, verbose=verbose)
	invisible(pred)
}

e0.predict.extra <- function(sim.dir=file.path(getwd(), 'bayesLife.output'), 
					prediction.dir=sim.dir, 
					countries = NULL, save.as.ascii=1000, verbose=TRUE, ...) {
	# Run prediction for given countries/regions (as codes). If they are not given it will be set to countries 
	# for which there are MCMC results but no prediction.
	# It is to be used after running run.e0.mcmc.extra
	
	mcmc.set <- get.e0.mcmc(sim.dir)
	if(is.null(mcmc.set))
		stop('Error in "sim.dir" argument.')
	pred <- get.e0.prediction(sim.dir=prediction.dir)
	if(is.null(pred))
		stop('Error in "prediction.dir" argument.')
	if(length(setdiff(pred$mcmc.set$meta$regions$country_code, mcmc.set$meta$regions$country_code)) > 0)
		stop('Prediction is inconsistent with the mcmc results. Use e0.predict.')
	if(is.null(countries)) {
		countries.idx <- (1:mcmc.set$meta$nr.countries)[!is.element(mcmc.set$meta$regions$country_code, 
												pred$mcmc.set$meta$regions$country_code)]
	} else {
		countries.idx <- (1:mcmc.set$meta$nr.countries)[is.element(mcmc.set$meta$regions$country_code,
												countries)]
	}
	if(length(countries.idx) == 0) {
		cat('\nNothing to be done.\n')
		return(invisible(pred))	
	}
	new.pred <- make.e0.prediction(mcmc.set, start.year=pred$start.year, end.year=pred$end.year, replace.output=FALSE,
									nr.traj=pred$nr.traj, burnin=pred$burnin,
									countries=countries.idx, save.as.ascii=0, output.dir=prediction.dir,
									force.creating.thinned.mcmc=TRUE,
									write.summary.files=FALSE, verbose=verbose)
									
	# merge the two predictions
	code.other.countries <- setdiff(pred$mcmc.set$meta$regions$country_code, 
									mcmc.set$meta$regions$country_code[countries.idx])
	idx.pred.others <- (1:pred$mcmc.set$meta$nr.countries)[is.element(pred$mcmc.set$meta$regions$country_code, 
												code.other.countries)]
	idx.other.countries <- (1:mcmc.set$meta$nr.countries)[is.element(mcmc.set$meta$regions$country_code,
												code.other.countries)]										
	prev.pred <- pred
	pred$quantiles <- new.pred$quantiles
	pred$quantiles[idx.other.countries,,] <- prev.pred$quantiles[idx.pred.others,,]
	
	pred$traj.mean.sd <- new.pred$traj.mean.sd
	pred$traj.mean.sd[idx.other.countries,,] <- prev.pred$traj.mean.sd[idx.pred.others,,]
	
	pred$e0.matrix.reconstructed <- new.pred$e0.matrix.reconstructed
	pred$e0.matrix.reconstructed[,idx.other.countries] <- prev.pred$e0.matrix.reconstructed[,idx.pred.others]
	
	pred$mcmc.set <- new.pred$mcmc.set
	
	if (has.e0.jmale.prediction(pred)) {
		bayesLife.prediction <- .do.e0.jmale.predict.extra(pred, countries.idx, idx.other.countries, idx.pred.others, ..., verbose=verbose)
		bayesTFR:::do.convert.trajectories(pred=get.e0.jmale.prediction(bayesLife.prediction), n=save.as.ascii, 
										output.dir=bayesLife.prediction$joint.male$output.directory, verbose=verbose)
	} else {
		# save updated prediction
		bayesLife.prediction <- pred
		prediction.file <- file.path(pred$output.directory, 'prediction.rda')
		save(bayesLife.prediction, file=prediction.file)
	}
	# convert trajectories and create summary files
	bayesTFR:::do.convert.trajectories(pred=bayesLife.prediction, n=save.as.ascii, output.dir=pred$output.directory, 
							verbose=verbose)
	do.write.e0.projection.summary(bayesLife.prediction, output.dir=pred$output.directory)
	
	cat('\nPrediction stored into', pred$output.directory, '\n')
	invisible(bayesLife.prediction)
}


make.e0.prediction <- function(mcmc.set, start.year=NULL, end.year=2100, replace.output=FALSE,
								nr.traj = NULL, thin=NULL, burnin=0, countries = NULL,
							    save.as.ascii=1000, output.dir = NULL, write.summary.files=TRUE, 
							    force.creating.thinned.mcmc=FALSE,
							    verbose=verbose){
	# if 'countries' is given, it is an index
	present.year <- if(is.null(start.year)) mcmc.set$meta$present.year else start.year - 5
	nr_project <- length(seq(present.year+5, end.year, by=5))
	cat('\nPrediction from', present.year+5, 'until', end.year, '(i.e.', nr_project, 'projections)\n\n')
			
	total.iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
	stored.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burnin)
	mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
	if(!is.null(nr.traj) && !is.null(thin)) {
		warning('Both nr.traj and thin are given. Argument thin will be ignored.')
		thin <- NULL
	}
	if(is.null(nr.traj)) nr.traj <- min(stored.iter, 2000)
	else {
		if (nr.traj > stored.iter) 
			warning('nr.traj is larger than the available MCMC sample. Only ', stored.iter, ' trajectories will be generated.')
		nr.traj <- min(nr.traj, stored.iter)	
	}
	if(is.null(thin)) thin <- floor(stored.iter/nr.traj * mcthin)
	if(stored.iter <= 0 || thin == 0)
		stop('The number of simulations is 0. Burnin might be larger than the number of simulated values, or # trajectories is too big.')

	#setup output directory
	if (is.null(output.dir)) output.dir <- mcmc.set$meta$output.dir
	outdir <- file.path(output.dir, 'predictions')

	if(is.null(countries)) {
		if(!replace.output && has.e0.prediction(sim.dir=output.dir))
			stop('Prediction in ', outdir,
				' already exists.\nSet replace.output=TRUE if you want to overwrite existing projections.')
		unlink(outdir, recursive=TRUE)
		write.to.disk <- TRUE
		if(!file.exists(outdir)) 
			dir.create(outdir, recursive=TRUE)
	} else write.to.disk <- FALSE
	
	thinned.mcmc <- get.thinned.e0.mcmc(mcmc.set, thin=thin, burnin=burnin)
	has.thinned.mcmc <- (!is.null(thinned.mcmc) && thinned.mcmc$meta$parent.iter == total.iter 
							&& mcmc.set$meta$nr.countries == thinned.mcmc$meta$nr.countries)
	unblock.gtk('bDem.e0pred')
	load.mcmc.set <- if(has.thinned.mcmc && !force.creating.thinned.mcmc) thinned.mcmc
		 			 else create.thinned.e0.mcmc(mcmc.set, thin=thin, burnin=burnin, 
					 					output.dir=output.dir, verbose=verbose)
	nr_simu <- load.mcmc.set$mcmc.list[[1]]$finished.iter	
	if (verbose) cat('Load world-level parameters.\n')
	var.par.names <- c('omega')
	# load only the first par to check if everything is o.k.
	var.par.values <- get.e0.parameter.traces(load.mcmc.set$mcmc.list, var.par.names, burnin=0)

	prediction.countries <- if(is.null(countries)) 1:mcmc.set$meta$nr.countries else countries
	nr_countries <- mcmc.set$meta$nr.countries
	e0.matrix.reconstructed <- get.e0.reconstructed(mcmc.set$meta$e0.matrix, mcmc.set$meta)
	#le0.matrix <- dim(e0.matrix.reconstructed)[1]
	present.year.index <- bayesTFR:::get.estimation.year.index(mcmc.set$meta, present.year)
	le0.matrix <- present.year.index
	
	quantiles.to.keep <- c(0,0.025,0.05,0.1,0.2,0.25,0.3,0.4,0.5,0.6,0.7,0.75,0.8,0.9,0.95,0.975,1)
	PIs_cqp <- array(NA, c(nr_countries, length(quantiles.to.keep), nr_project+1))
	dimnames(PIs_cqp)[[2]] <- quantiles.to.keep
	proj.middleyears <- bayesTFR:::get.prediction.years(mcmc.set$meta, nr_project+1, present.year.index)
	dimnames(PIs_cqp)[[3]] <- proj.middleyears
	mean_sd <- array(NA, c(nr_countries, 2, nr_project+1))

	var.par.names.cs <- c('Triangle.c', 'k.c', 'z.c')
	
	country.counter <- 0
	status.for.gui <- paste('out of', length(prediction.countries), 'countries.')
	gui.options <- list()
	#########################################
	for (country in prediction.countries){
	#for (country in c(23)){
	#########################################
		if(getOption('bDem.e0pred', default=FALSE)) {
			# This is to unblock the GUI, if the run is invoked from bayesDem
			# and pass info about its status
			# In such a case the gtk libraries are already loaded
			country.counter <- country.counter + 1
			gui.options$bDem.e0pred.status <- paste('finished', country.counter, status.for.gui)
			unblock.gtk('bDem.e0pred', gui.options)
		}

		country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
		if (verbose) {			
 			cat('e0 projection for country', country, country.obj$name, 
 						'(code', country.obj$code, ')\n')
 		}
 		if (is.element(country.obj$code, load.mcmc.set$meta$regions$country_code)) {
			cs.par.values <- get.e0.parameter.traces.cs(load.mcmc.set$mcmc.list, country.obj, 
								var.par.names.cs, burnin=0)
		} else { # there are no thinned traces for this country, use the full traces
			cs.par.values <- get.e0.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, 
								var.par.names.cs, burnin=burnin)
			selected.simu <- bayesTFR:::get.thinning.index(nr_simu, dim(cs.par.values)[1])
			if (length(selected.simu$index) < nr_simu)
				selected.simu$index <- sample(selected.simu$index, nr_simu, replace=TRUE)
			cs.par.values <- cs.par.values[selected.simu$index,]
		}
		all.e0 <- mcmc.set$meta$e0.matrix[, country]
		lall.e0 <- length(all.e0)
		this.Tc_end <-  mcmc.set$meta$Tc.index[[country]][length(mcmc.set$meta$Tc.index[[country]])]
		this.T_end <- min(this.Tc_end, le0.matrix)
		nmissing <- le0.matrix - this.T_end
		missing <- (this.T_end+1):le0.matrix

		if (verbose && nmissing > 0) 
			cat('\t', nmissing, 'data points reconstructed.\n')

		this.nr_project <- nr_project + nmissing
		#sum.delta <- apply(cs.par.values[,1:4], 1, sum)
		#use.traj <- which(sum.delta <= 110)
		trajectories <- matrix(NA, this.nr_project+1, nr_simu)
		#trajectories <- matrix(NA, this.nr_project+1, length(use.traj))
    	for(j in 1:nr_simu){
    	#for(j in 1:length(use.traj)){
    		#k <- use.traj[j]
    		trajectories[1,j] <- all.e0[this.T_end]
    		if(nmissing == 0 && this.Tc_end > le0.matrix) { # use observed data on projection spots
    			#stop('')
    			trajectories[2:(lall.e0 - le0.matrix+1),j]<- all.e0[(le0.matrix+1):lall.e0]
    			proj.idx <- (lall.e0 - le0.matrix + 2):(this.nr_project+1)
    			last.val.idx <- lall.e0
    		} else {
    			proj.idx <- 2:(this.nr_project+1)
    			last.val.idx <- this.T_end
    		}
           trajectories[proj.idx,j]<-e0.proj.le.SDPropToLoess(cs.par.values[j,], 
           							all.e0[last.val.idx], 
           							kap=var.par.values[j,'omega'],n.proj=length(proj.idx),
           							p1=mcmc.set$meta$dl.p1, p2=mcmc.set$meta$dl.p2, 
           							const.var=mcmc.set$meta$constant.variance)
    	}
    	if (nmissing > 0) {
    		e0.matrix.reconstructed[(this.T_end+1):le0.matrix,country] <- apply(matrix(trajectories[2:(nmissing+1),],
    											 nrow=nmissing), 1, quantile, 0.5, na.rm = TRUE)
    		trajectories <- trajectories[(nmissing+1):nrow(trajectories),]
    		trajectories[1,] <- quantile(trajectories[1,], 0.5, na.rm = TRUE)
    	}
    	#stop('')
		save(trajectories, file = file.path(outdir, paste('traj_country', country.obj$code, '.rda', sep='')))
 		PIs_cqp[country,,] = apply(trajectories, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		mean_sd[country,1,] <- apply(trajectories, 1, mean, na.rm = TRUE)
 		mean_sd[country,2,] = apply(trajectories, 1, sd, na.rm = TRUE)
	}
	mcmc.set <- remove.e0.traces(mcmc.set)
	bayesLife.prediction <- structure(list(
				quantiles = PIs_cqp,
				traj.mean.sd = mean_sd,
				nr.traj=nr_simu,
				e0.matrix.reconstructed = e0.matrix.reconstructed,
				output.directory=outdir,
				mcmc.set=load.mcmc.set,
				nr.projections=nr_project,
				burnin=burnin,
				end.year=end.year,
				start.year=start.year,
				present.year.index=present.year.index,
				present.year.index.all=present.year.index + (
						if(!is.null(mcmc.set$meta$suppl.data$regions)) nrow(mcmc.set$meta$suppl.data$e0.matrix) else 0)
				),
				class='bayesLife.prediction')
		
	if(write.to.disk) {		
		prediction.file <- file.path(outdir, 'prediction.rda')
		save(bayesLife.prediction, file=prediction.file)
	
		bayesTFR:::do.convert.trajectories(pred=bayesLife.prediction, n=save.as.ascii, output.dir=outdir, 
										verbose=verbose)
		if(write.summary.files)
			do.write.e0.projection.summary(bayesLife.prediction, output.dir=outdir)
	
		cat('\nPrediction stored into', outdir, '\n')
	}
	invisible(bayesLife.prediction)
}

remove.e0.traces <- function(mcmc.set) {
	for (i in 1:length(mcmc.set$mcmc.list)) 
		mcmc.set$mcmc.list[[i]]$traces <- 0
	invisible(mcmc.set)
}

get.projection.summary.header.bayesLife.prediction <- function(pred, ...) 
		return (list(revision='RevID', variant='VarID', country='LocID', year='TimeID', indicator='IndicatorID', sex='SexID', tfr='Value'))
		
get.UN.variant.names.bayesLife.prediction <- function(pred, ...) 
		return(c('BHM median', 'BHM80 lower',  'BHM80 upper', 'BHM95 lower',  'BHM95 upper', 'Constant mortality'))
	
get.friendly.variant.names.bayesLife.prediction <- function(pred, ...)
	return(c('median', 'lower 80', 'upper 80', 'lower 95', 'upper 95', 'constant'))	

convert.e0.trajectories <- function(dir=file.path(getwd(), 'bayesLife.output'), 
								 n=1000, output.dir=NULL, 
								 verbose=FALSE) {
	# Converts all trajectory rda files into UN ascii, selecting n trajectories by equal spacing.
	if(n <= 0) return()
	pred <- get.e0.prediction(sim.dir=dir)
	predsex <- pred$mcmc.set$meta$sex
	preds <- list()
	preds[[predsex]] <- pred
	if(has.e0.jmale.prediction(pred)) preds[['M']] <- get.e0.jmale.prediction(pred)
	for (sex in c('F', 'M')) {
		if(is.null(preds[[sex]])) next
		if (is.null(output.dir)) outdir <- preds[[sex]]$output.directory
		else {
			if(length(preds) > 1) outdir <- file.path(output.dir,sex)
			else outdir <- output.dir
		}
		if(!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
		cat('Converting ', list(M='Male', F='Female')[[sex]], ' trajectories from', dir, '\n')
		bayesTFR:::do.convert.trajectories(pred=preds[[sex]], n=n, output.dir=outdir, verbose=verbose)
	}
}

write.e0.projection.summary <- function(dir=file.path(getwd(), 'bayesLife.output'), 
									 output.dir=NULL, revision=NULL, adjusted=FALSE) {
# Writes four prediction summary files, one in a user-friendly format, one in a UN-format, one for each sex.
	pred <- get.e0.prediction(sim.dir=dir)
	predsex <- pred$mcmc.set$meta$sex
	preds <- list()
	preds[[predsex]] <- pred
	if(has.e0.jmale.prediction(pred)) preds[['M']] <- get.e0.jmale.prediction(pred)
	for (sex in c('F', 'M')) {
		if(is.null(preds[[sex]])) next
		if (is.null(output.dir)) outdir <- preds[[sex]]$output.directory
		else {
			if(length(preds) > 1) outdir <- file.path(output.dir,sex)
			else outdir <- output.dir
		}
		if(!file.exists(outdir)) dir.create(outdir, recursive=TRUE)
		do.write.e0.projection.summary(preds[[sex]], outdir, sex=sex, revision=revision, adjusted=adjusted)
	}
}
		
do.write.e0.projection.summary <- function(pred, output.dir, sex=NULL, revision=NULL, adjusted=FALSE) {
	if (is.null(sex)) sex <- pred$mcmc.set$meta$sex
	bayesTFR:::do.write.projection.summary(pred, output.dir, revision=revision, indicator.id=10, 
				sex.id=c(M=1,F=2)[sex], adjusted=adjusted)
}
				
get.traj.ascii.header.bayesLife.mcmc.meta <- function(meta, ...) 
	return (list(country_code='LocID', period='Period', year='Year', trajectory='Trajectory', tfr='e0'))
	
get.data.imputed.bayesLife.prediction <- function(pred, ...)
	return(get.e0.reconstructed(pred$e0.matrix.reconstructed, pred$mcmc.set$meta))
	
get.data.imputed.for.country.bayesLife.prediction <- function(pred, country.index, ...)
	return(bayesTFR:::get.observed.with.supplemental(country.index, pred$e0.matrix.reconstructed, 
					pred$mcmc.set$meta$suppl.data, 'e0.matrix'))
	
get.e0.reconstructed <- function(data, meta) {
	return(if(is.null(data)) meta$e0.matrix.all else data)
}

e0.median.reset <- function(sim.dir, countries, joint.male=FALSE) {
	for(country in countries) pred <- e0.median.shift(sim.dir, country, reset=TRUE, joint.male=joint.male)
	invisible(pred)
}

get.e0.shift <- function(country.code, pred) return(bayesTFR::get.tfr.shift(country.code, pred))

e0.median.shift <- function(sim.dir, country, reset=FALSE, shift=0, from=NULL, to=NULL, joint.male=FALSE) {
	pred <- get.e0.prediction(sim.dir, joint.male=joint.male)
	new.pred <- bayesTFR:::.bdem.median.shift(pred, type='e0', country=country, reset=reset, 
				shift=shift, from=from, to=to)
	if(joint.male) {
		predF <- get.e0.prediction(sim.dir)
		predF$joint.male <- new.pred
		store.bayesLife.prediction(predF)
	} else store.bayesLife.prediction(new.pred)
	invisible(new.pred)
}

e0.median.set <- function(sim.dir, country, values, years=NULL, joint.male=FALSE) {
	pred <- get.e0.prediction(sim.dir, joint.male=joint.male)
	new.pred <- bayesTFR:::.bdem.median.set(pred, type='e0', country=country, 
								values=values, years=years)
	if(joint.male) {
		predF <- get.e0.prediction(sim.dir)
		predF$joint.male <- new.pred
		store.bayesLife.prediction(predF)
	} else store.bayesLife.prediction(new.pred)
	invisible(new.pred)
}

e0.jmale.estimate <- function(mcmc.set, countries.index=NULL, 
								estDof.eq1 = TRUE, start.eq1 = list(dof = 2), 
								max.e0.eq1 = 83, estDof.eq2 = TRUE, start.eq2 = list(dof = 2), 
								constant.gap.eq2=TRUE, my.e0.file=NULL, my.locations.file=NULL, verbose=FALSE) {
	# Estimate coefficients for joint prediction of female and male e0
	unblock.gtk('bDem.e0pred', list(bDem.e0pred.status='estimating joint male'))
	if (is.null(countries.index)) countries.index <- 1:get.nr.countries.est(mcmc.set$meta)
	e0f.data <- get.data.matrix(mcmc.set$meta)[,countries.index]
	e0m.data <- get.wpp.e0.data.for.countries(mcmc.set$meta, sex='M', my.e0.file=my.e0.file,
					my.locations.file=my.locations.file, verbose=verbose)$e0.matrix[,countries.index]
	T <- dim(e0f.data)[1] - 1 
	if(verbose) {
		cat('\nEstimating coefficients for joint female and male prediction.')
		cat('\nUsing', length(countries.index), 'countries and', T+1, 'time periods.\n\n')
	}
	G <- e0f.data - e0m.data # observed gap

	first.year <- max(1953, as.integer(rownames(e0f.data)[1])) # in case start.year is later than 1953
	if(first.year > 1953)
		warning("Data for 1950-1955 not available. Estimation of the gap model may not be correct.")
	e0.1953 <- rep(e0f.data[as.character(first.year),], each=T)
	dep.var <- as.numeric(G[2:nrow(G),])
	e0F <- as.numeric(e0f.data[2:(T+1),])
	data <- data.frame(
					G=dep.var, # dependent variable
					# covariates
					e0.1953=e0.1953, 
					Gprev=as.numeric(G[1:T,]),
					e0=e0F,
					e0d75=pmax(0,e0F-75)
	)
	use.eq1 <- !is.na(dep.var) & !is.na(e0F) & e0F <= max.e0.eq1
	data.eq1 <- data[use.eq1,]
	fit.eq1 <- tlm(G~., data=data.eq1, estDof = estDof.eq1, start=start.eq1)
	if(verbose)
		print(summary(fit.eq1))
	errscale.eq1<-as.numeric(exp(fit.eq1$scale.fit$coefficients[1]))
	errsd.eq1<-sqrt(errscale.eq1)

	use.eq2 <- !is.na(dep.var) & !is.na(e0F) & e0F > max.e0.eq1
	data.eq2 <- data[use.eq2,]
	if(verbose) 
		cat('\n\nUsing', nrow(data.eq2), 'data points for equation 2.\n\n')
	if(!constant.gap.eq2) {
		fit.eq2 <- tlm(G~-1+Gprev, data=data.eq2, start = start.eq2, estDof = estDof.eq2)
		if(verbose)
			print(summary(fit.eq2))
		errscale.eq2<-as.numeric(exp(fit.eq2$scale.fit$coefficients[1]))
		errsd.eq2<-sqrt(errscale.eq2)
		coef2 <- fit.eq2$loc.fit$coefficients
		dof2 <- fit.eq2$dof
	} else {# constant gap in eq.2
		dof2 <- NULL
		coef2 <- c(Gprev=1)
		errsd.eq2 <- sqrt(mean((data.eq2$G - data.eq2$Gprev)^2))
	}
	return(list(eq1 = list(coefficients=fit.eq1$loc.fit$coefficients, 
						   sigma=errsd.eq1, dof = fit.eq1$dof),
				eq2 = list(coefficients=coef2, sigma=errsd.eq2, dof = dof2)
				))
}

e0.jmale.predict <- function(e0.pred, estimates=NULL, gap.lim=c(0,18),  #gap.lim.eq2=c(3,9),	
								max.e0.eq1.pred=83, my.e0.file=NULL, my.locations.file=NULL, 
								save.as.ascii=1000, verbose=TRUE, ...) {
	# Predicting male e0 from female predictions. estimates is the result of 
	# the e0.jmale.estimate function. If it is NULL, the estimation is performed 
	# using the ... arguments
	# If my.e0.file given, it should be a male e0 file. 
	
	meta <- e0.pred$mcmc.set$meta
	if (meta$sex != 'F') stop('The prediction object must be a result of FEMALE projections.')
	if(is.null(estimates)) 
		estimates <- e0.jmale.estimate(e0.pred$mcmc.set, verbose=verbose, 
								my.e0.file=my.e0.file, my.locations.file=my.locations.file, ...)

	e0mwpp <- get.wpp.e0.data.for.countries(meta, sex='M', my.e0.file=my.e0.file, 
											my.locations.file=my.locations.file, verbose=verbose)
	e0m.data <- e0mwpp$e0.matrix
	meta.changes <- list(sex='M', e0.matrix=e0m.data, e0.matrix.all=e0mwpp$e0.matrix.all, suppl.data=e0mwpp$suppl.data)
	meta.changes$Tc.index <- .get.Tcindex(meta.changes$e0.matrix, cnames=meta$regions$country_name)
	if(!is.null(meta.changes$suppl.data$e0.matrix))
		meta.changes$suppl.data$Tc.index <- .get.Tcindex(meta.changes$suppl.data$e0.matrix, stop.if.less.than2=FALSE)
	prediction.file <- file.path(e0.pred$output.directory, 'prediction.rda')
	joint.male <- e0.pred
	joint.male$output.directory <- file.path(e0.pred$output.directory, 'joint_male')
	joint.male$e0.matrix.reconstructed <- e0m.data
	joint.male$fit <- estimates
	joint.male$meta.changes <- meta.changes
	joint.male$mcmc.set <- NULL
	joint.male$joint.male <- NULL
	joint.male$pred.pars <- list(gap.lim=gap.lim, #gap.lim.eq2=gap.lim.eq2, 
								max.e0.eq1.pred=max.e0.eq1.pred)
	
	if(file.exists(joint.male$output.directory)) unlink(joint.male$output.directory, recursive=TRUE)
	dir.create(joint.male$output.directory, recursive=TRUE)
	bayesLife.prediction <- .do.jmale.predict(e0.pred, joint.male, 1:get.nr.countries(meta),  
								gap.lim=gap.lim, #gap.lim.eq2=gap.lim.eq2, 
								eq2.age.start=max.e0.eq1.pred, verbose=verbose)
	save(bayesLife.prediction, file=prediction.file)
	cat('\nPrediction stored into', joint.male$output.directory, '\n')
	bayesTFR:::do.convert.trajectories(pred=get.e0.jmale.prediction(bayesLife.prediction), n=save.as.ascii, 
										output.dir=joint.male$output.directory, verbose=verbose)
	do.write.e0.projection.summary(get.e0.jmale.prediction(bayesLife.prediction), output.dir=joint.male$output.directory, sex='M')
	invisible(bayesLife.prediction)
}

.do.e0.jmale.predict.extra <- function(e0.pred, countries.idx, idx.other.to.new, idx.other.to.old,
									gap.lim.eq1=c(0,18),  gap.lim.eq2=c(3,9), max.e0.eq1.pred=83, my.e0.file=NULL, 
									my.locations.file=NULL, verbose=TRUE) {
	# called from e0.predict.extra
	if (!has.e0.jmale.prediction(e0.pred)) stop('Joint female-male prediction must be available for e0.pred. Use e0.jmale.predict.')
	joint.male <- get.e0.jmale.prediction(e0.pred)
	meta <- e0.pred$mcmc.set$meta # from female sim
	male.meta <- joint.male$meta.changes
	meta$sex <- "M"
	e0mwpp <- set.e0.wpp.extra(meta, meta$regions$country_code[countries.idx], my.e0.file=my.e0.file, my.locations.file=my.locations.file)
	# merge e0 matrices
	for(mat in c("e0.matrix", "e0.matrix.all")) {
		tmp <- meta[[mat]] # the right size
		tmp[,idx.other.to.new] <- male.meta[[mat]][,idx.other.to.old] # replace by the male original (non-extra) data
		tmp[,countries.idx] <- e0mwpp[[mat]] # replace by data from extra countries
		joint.male$meta.changes[[mat]] <- tmp
	}
	#e0mwpp <- get.wpp.e0.data.for.countries(meta, sex='M', my.e0.file=my.e0.file, verbose=verbose)
	#e0m.data <- e0mwpp$e0.matrix
	#e0m.data[,idx.other.to.new] <- joint.male$meta.changes$e0.matrix[,idx.other.to.old]
	#e0mwpp$e0.matrix.all[,idx.other.to.new] <- joint.male$meta.changes$e0.matrix.all[,idx.other.to.old]
	#meta.changes <- list(sex='M', e0.matrix=e0m.data, e0.matrix.all=e0mwpp$e0.matrix.all, suppl.data=e0mwpp$suppl.data)
	Tc.index <- .get.Tcindex(joint.male$meta.changes$e0.matrix, cnames=meta$regions$country_name)
	#meta.changes$Tc.index <- joint.male$meta.changes$Tc.index 
	joint.male$meta.changes$Tc.index[countries.idx] <- Tc.index[countries.idx]
	# We don't need to get the supplemental data for extra countries since they are not used for estimation
	#Tc.index <- .get.Tcindex(meta.changes$suppl.data$e0.matrix, stop.if.less.than2=FALSE)
	#meta.changes$suppl.data$Tc.index <- joint.male$meta.changes$suppl.data$Tc.index
	#meta.changes$suppl.data$Tc.index[countries.idx] <- Tc.index[countries.idx]
	
	prediction.file <- file.path(e0.pred$output.directory, 'prediction.rda')
	#joint.male$meta.changes <- meta.changes
	reconstructed <- joint.male$meta.changes$e0.matrix
	reconstructed[1:nrow(joint.male$e0.matrix.reconstructed),1:ncol(joint.male$e0.matrix.reconstructed)] <- joint.male$e0.matrix.reconstructed
	joint.male$e0.matrix.reconstructed <- reconstructed
	new.pred <- .do.jmale.predict(e0.pred, joint.male, countries.idx, gap.lim=joint.male$pred.pars$gap.lim, 
									#gap.lim.eq2=joint.male$gap.lim.eq2,
									eq2.age.start=joint.male$pred.pars$max.e0.eq1.pred, verbose=verbose)
	new.jmale <- new.pred$joint.male
	prev.jmale <- e0.pred$joint.male
	joint.male$quantiles <- new.jmale$quantiles
	joint.male$quantiles[idx.other.to.new,,] <- prev.jmale$quantiles[idx.other.to.old,,]
	
	joint.male$traj.mean.sd <- new.jmale$traj.mean.sd
	joint.male$traj.mean.sd[idx.other.to.new,,] <- prev.jmale$traj.mean.sd[idx.other.to.old,,]
	
	joint.male$e0.matrix.reconstructed <- new.jmale$e0.matrix.reconstructed
	joint.male$e0.matrix.reconstructed[,idx.other.to.new] <- prev.jmale$e0.matrix.reconstructed[,idx.other.to.old]

	
	bayesLife.prediction <- e0.pred
	bayesLife.prediction$joint.male <- joint.male
	save(bayesLife.prediction, file=file.path(e0.pred$output.directory, 'prediction.rda'))
	cat('\nPrediction stored into', joint.male$output.directory, '\n')
	invisible(bayesLife.prediction)
}


.do.jmale.predict <- function(e0.pred, joint.male, countries, gap.lim, #gap.lim.eq2, 
								eq2.age.start=NULL, verbose=FALSE) {
	predict.one.trajectory <- function(Gprev, ftraj) {
		mtraj <- rep(NA, length(ftraj))						
		for(time in 1:length(ftraj)) {
			if(ftraj[time] <= maxe0) { # 1st part of Equation 3.1
				Gtdeterm <- (estimates$eq1$coefficients[1] + # intercept
				   			 estimates$eq1$coefficients['Gprev']*Gprev +
				   			 estimates$eq1$coefficients['e0.1953']*e0f.data[first.year,icountry] +
				   			 estimates$eq1$coefficients['e0']*ftraj[time] +
					   		 estimates$eq1$coefficients['e0d75']*max(0, ftraj[time]-75))
				Gt <- Gtdeterm + estimates$eq1$sigma*rt(1,estimates$eq1$dof)
				while(Gt < gap.lim[1] || Gt > gap.lim[2]) 
					Gt <- Gtdeterm + estimates$eq1$sigma*rt(1,estimates$eq1$dof)
			} else {  # 2nd part of Equation 3.1
				Gtdeterm <- estimates$eq2$coefficients['Gprev']*Gprev
				error <- if(is.null(estimates$eq2$dof)) rnorm(1, sd=estimates$eq2$sigma) 
			    			else estimates$eq2$sigma*rt(1,estimates$eq2$dof)
				Gt <- Gtdeterm + error					
				while(Gt < gap.lim[1] || Gt > gap.lim[2]) {
					Gt <- Gtdeterm + if(is.null(estimates$eq2$dof)) rnorm(1, sd=estimates$eq2$sigma) 
							else estimates$eq2$sigma*rt(1,estimates$eq2$dof)
				}
			}
			mtraj[time] <- ftraj[time] - Gt
			Gprev <- Gt
		}
		return(mtraj)
	}							
									
	unblock.gtk('bDem.e0pred', list(bDem.e0pred.status='predicting joint male'))
	bayesLife.prediction <- e0.pred
	bayesLife.prediction$joint.male <- joint.male
	meta <- e0.pred$mcmc.set$meta
	quantiles <- array(NA, dim(e0.pred$quantiles))
	dimnames(quantiles) <- dimnames(e0.pred$quantiles)
	traj.mean.sd <- array(NA, dim(e0.pred$traj.mean.sd))
	dimnames(traj.mean.sd) <- dimnames(e0.pred$traj.mean.sd)
	e0f.data <- get.e0.reconstructed(e0.pred$e0.matrix.reconstructed, meta)
	maxe0 <- if(is.null(eq2.age.start)) max(e0f.data) else eq2.age.start
	e0m.data <- joint.male$e0.matrix.reconstructed
	quantiles.to.keep <- as.numeric(dimnames(e0.pred$quantiles)[[2]])
	first.year <- max(1953, as.integer(rownames(e0f.data)[1]))
	if(first.year > 1953)
		warning("Data for 1950-1955 not available. Projection of the gap model may not be correct.")
	first.year <- as.character(first.year)
	estimates <- joint.male$fit
	for (icountry in countries) {
		country <- get.country.object(icountry, meta, index=TRUE)
		if(verbose)
			cat('\ne0 male projection for country', icountry, country$name, 
 						'(code', country$code, ')')
		trajectoriesF <- bayesTFR:::get.trajectories(e0.pred, country$code)$trajectories
		Mtraj <- matrix(NA, nrow=nrow(trajectoriesF), ncol=ncol(trajectoriesF))
		#G1 <- e0f.data[Tc[icountry],icountry] - e0m.data[Tc[icountry],icountry]
		Tc <- joint.male$meta.changes$Tc.index[[icountry]][length(joint.male$meta.changes$Tc.index[[icountry]])]
		G1 <- e0f.data[Tc,icountry] - e0m.data[Tc,icountry]
		last.obs.index <- if(is.null(e0.pred$present.year.index)) nrow(e0m.data) else e0.pred$present.year.index
		if(Tc < last.obs.index) { # imputing data
			imp.index <- (last.obs.index - Tc + 1):last.obs.index
			e0m.data[imp.index,icountry] <- predict.one.trajectory(G1, e0f.data[imp.index, icountry])
			G1 <- e0f.data[last.obs.index,icountry] - e0m.data[last.obs.index,icountry]
		}
		for (itraj in 1:dim(trajectoriesF)[2]) {
			Mtraj[1,itraj] <- e0m.data[last.obs.index,icountry]
			Mtraj[2:nrow(Mtraj),itraj] <- predict.one.trajectory(G1, trajectoriesF[2:nrow(trajectoriesF),itraj])
		}
		quantiles[icountry,,] = apply(Mtraj, 1, quantile, quantiles.to.keep, na.rm = TRUE)
 		traj.mean.sd[icountry,1,] <- apply(Mtraj, 1, mean, na.rm = TRUE)
 		traj.mean.sd[icountry,2,] = apply(Mtraj, 1, sd, na.rm = TRUE)
 		trajectories <- Mtraj
		save(trajectories, file = file.path(joint.male$output.directory, 
								paste('traj_country', country$code, '.rda', sep='')))
	}
	bayesLife.prediction$joint.male$quantiles <- quantiles
	bayesLife.prediction$joint.male$traj.mean.sd <- traj.mean.sd
	bayesLife.prediction$joint.male$e0.matrix.reconstructed <- e0m.data
	return(bayesLife.prediction)
	
	
}

get.e0.jmale.prediction <- function(e0.pred) {
	male.pred <- e0.pred$joint.male
	if(is.null(male.pred)) stop('A joint male prediction does not exist for the given object. Use e0.jmale.predict to simulate male projections from existing female projections.')
	male.pred$mcmc.set <- e0.pred$mcmc.set
	for(item in names(male.pred$meta.changes))
		male.pred$mcmc.set$meta[[item]] <- male.pred$meta.changes[[item]]
	return(male.pred)
}

has.e0.jmale.prediction <- function(e0.pred) return(!is.null(e0.pred$joint.male))
