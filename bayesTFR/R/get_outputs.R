get.tfr.mcmc <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), chain.ids=NULL, low.memory=TRUE,
							burnin=0, verbose=FALSE) {
	############
	# Returns an object of class bayesTFR.mcmc.set
	############
	mcmc.file.path <- file.path(sim.dir, 'bayesTFR.mcmc.meta.rda')
	if(!file.exists(mcmc.file.path)) {
		warning('File ', mcmc.file.path, ' does not exist.')
		return(NULL)
	}
	load(file=mcmc.file.path)
	bayesTFR.mcmc.meta$output.dir <- normalizePath(sim.dir)
	if (is.null(chain.ids)) {
		mc.dirs.short <- list.files(sim.dir, pattern='^mc[0-9]+', full.names=FALSE)
		chain.ids <- as.integer(substring(mc.dirs.short, 3))
	} else {
		mc.dirs.short <- paste('mc', chain.ids, sep='')
	}
	ord.idx <- order(chain.ids)
	mc.dirs.short <- mc.dirs.short[ord.idx]
	chain.ids <- chain.ids[ord.idx]
	mcmc.chains <- list()
	counter<-1
	for (imc.d in chain.ids) {
		if (verbose)
			cat('Loading chain', imc.d, 'from disk. ')
		bayesTFR.mcmc <- local({
			load(file=file.path(sim.dir, mc.dirs.short[counter], 'bayesTFR.mcmc.rda'))
			bayesTFR.mcmc})
		mc <- c(bayesTFR.mcmc, list(meta=bayesTFR.mcmc.meta))
		class(mc) <- class(bayesTFR.mcmc)
		if (!low.memory) { # load full mcmc traces
			th.burnin <- get.thinned.burnin(mc, burnin)
			mc$traces <- load.tfr.parameter.traces.all(mc, burnin=th.burnin)
			mc$traces.burnin <- th.burnin
		} else { # traces will be loaded as they are needed
			mc$traces <- 0
			mc$traces.burnin <- 0
		}
		mc$output.dir <- mc.dirs.short[counter]
		if (verbose)
			cat('(mcmc.list[[', counter, ']]).\n')
		mcmc.chains[[counter]] <- mc
		counter <- counter+1
	}
	names(mcmc.chains) <- chain.ids
	return(structure(list(meta=bayesTFR.mcmc.meta, 
				mcmc.list=mcmc.chains), class='bayesTFR.mcmc.set'))
}

get.tfr3.mcmc <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), ...) {
	parent.mc <- get.tfr.mcmc(sim.dir)
	mc <- get.tfr.mcmc(file.path(sim.dir, 'phaseIII'), ...)
	mc$meta$parent <- parent.mc$meta
	mc$meta$regions <- parent.mc$meta$regions
	if(length(mc$mcmc.list) > 0) {
		for(chain in 1:length(mc$mcmc.list)) {
			mc$mcmc.list[[chain]]$meta <- mc$meta
		}
	}
	return(mc)					
}

has.tfr.mcmc <- function(sim.dir) {
	return(file.exists(file.path(sim.dir, 'bayesTFR.mcmc.meta.rda')))
}

has.tfr3.mcmc <- function(sim.dir) {
	return(has.tfr.mcmc(file.path(sim.dir, 'phaseIII')))
}

tfr.mcmc <- function(mcmc.set, chain.id) return (mcmc.set$mcmc.list[[chain.id]])

tfr.mcmc.list <- function(mcmc.set, chain.ids=NULL) {
	if(is.null(chain.ids)) return(mcmc.set$mcmc.list)
	return(mcmc.set$mcmc.list[chain.ids])
}

get.thinned.tfr.mcmc <- function(mcmc.set, thin=1, burnin=0) {
	dir.name <- file.path(mcmc.set$meta$output.dir, paste('thinned_mcmc', thin, burnin, sep='_'))
	if(file.exists(dir.name)) return(get.tfr.mcmc(dir.name))
	return(NULL)
}
	
create.thinned.tfr.mcmc <- function(mcmc.set, thin=1, burnin=0, output.dir=NULL, verbose=TRUE) {
	#Return a thinned mcmc.set object with burnin removed and all chanins collapsed into one
	mcthin <- max(sapply(mcmc.set$mcmc.list, function(x) x$thin))
	thin <- max(c(thin, mcthin))
	meta <- mcmc.set$meta
	total.iter <- get.stored.mcmc.length(mcmc.set$mcmc.list, burnin=burnin)
	meta$is.thinned <- TRUE
	meta$parent.iter <- get.total.iterations(mcmc.set$mcmc.list, burnin)
	meta$parent.meta <- mcmc.set$meta
	meta$nr.chains <- 1
	
	if(verbose) cat('\nStoring thinned mcmc:')
	# store the meta object
	meta$output.dir <- file.path(
		if(is.null(output.dir)) meta$output.dir else output.dir, 
			paste('thinned_mcmc', thin, burnin, sep='_'))
	if(!file.exists(meta$output.dir)) 
		dir.create(meta$output.dir, recursive=TRUE)
	store.bayesTFR.meta.object(meta, meta$output.dir)
	
	thin.index <- if(thin > mcthin) unique(round(seq(1, total.iter, by=thin/mcthin))) else 1:total.iter
	nr.points <- length(thin.index)
	
	#create one collapsed mcmc
	thinned.mcmc <- mcmc.set$mcmc.list[[1]]
	thinned.mcmc$meta <- meta
	thinned.mcmc$thin <- 1
	thinned.mcmc$id <- 1
	thinned.mcmc$traces <- 0
	thinned.mcmc$length <- nr.points
	thinned.mcmc$finished.iter <- nr.points
	thinned.mcmc$compression.type <- meta$compression.type
	thinned.mcmc$output.dir <- 'mc1'	
	outdir.thin.mcmc <- file.path(meta$output.dir, 'mc1')
	if(!file.exists(outdir.thin.mcmc)) dir.create(outdir.thin.mcmc)

	store.bayesTFR.object(thinned.mcmc, outdir.thin.mcmc)
	
	if(verbose) cat('\nStoring country-independent parameters ...')
	for (par in tfr.parameter.names(trans=FALSE)) {
		values <- get.tfr.parameter.traces(mcmc.set$mcmc.list, par, burnin,
											thinning.index=thin.index)
		write.values.into.file.cindep(par, values, outdir.thin.mcmc, compression.type=thinned.mcmc$compression.type)
	}
	if(verbose) cat('done.\nStoring country-specific parameters ...')
	par.names.cs <- tfr.parameter.names.cs(trans=FALSE)
	for (country in mcmc.set$meta$id_DL){
		country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
		for (par in par.names.cs) {
			values <- get.tfr.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
											burnin=burnin, thinning.index=thin.index)
			write.values.into.file.cdep(par, values, outdir.thin.mcmc, country.code=country.obj$code,
										compression.type=thinned.mcmc$compression.type)
		}
	}
	if (mcmc.set$meta$nr_countries > mcmc.set$meta$nr_countries_estimation) {
		.update.thinned.extras(mcmc.set, (mcmc.set$meta$nr_countries_estimation+1):mcmc.set$meta$nr_countries,
								burnin=burnin, nr.points=nr.points, dir=outdir.thin.mcmc, verbose=verbose)
	}
	invisible(structure(list(meta=meta, mcmc.list=list(thinned.mcmc)), class='bayesTFR.mcmc.set'))
}

.update.thinned.extras <- function (mcmc.set, country.index, burnin, nr.points, dir, verbose=TRUE) {
	if(verbose) cat('done.\nStoring country-specific parameters for extra countries ...')
	# thin mcmc for extra countries (they can have different length than the other countries)
	par.names.cs <- tfr.parameter.names.cs(trans=FALSE)
	for (country in country.index){
		if(!(country %in% mcmc.set$meta$id_DL)) next
		country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
		for (par in par.names.cs) {
			values <- get.tfr.parameter.traces.cs(mcmc.set$mcmc.list, country.obj, par, 
											burnin=burnin)
			selected.simu <- get.thinning.index(nr.points, dim(values)[1])
			if (length(selected.simu$index) < nr.points)
				selected.simu$index <- sample(selected.simu$index, nr.points, replace=TRUE)
			values <- values[selected.simu$index,]
			write.values.into.file.cdep(par, values, dir, country.code=country.obj$code, 
								compression.type=mcmc.set$meta$compression.type)
		}
	}
	if(verbose) cat('done.\n')
}

.update.meta.for.thinned.mcmc <- function(thin.mcmc.set, new.mcmc.set) {
	# Updating meta object when extra countries are added
	keep.meta <- thin.mcmc.set$meta
	thin.meta <- new.mcmc.set$meta
	for(item in c('output.dir', 'is.thinned', 'parent.iter', 'nr.chains', 'compression.type'))
		thin.meta[[item]] <- keep.meta[[item]]
	thin.meta$parent.meta <- new.mcmc.set$meta
	for (ichain in 1:length(thin.mcmc.set$mcmc.list)) {
		thin.mcmc.set$mcmc.list[[ichain]]$meta <- thin.meta
		store.bayesTFR.object(thin.mcmc.set$mcmc.list[[ichain]], 
				file.path(thin.meta$output.dir, thin.mcmc.set$mcmc.list[[ichain]]$output.dir))
	}
	store.bayesTFR.meta.object(thin.meta, thin.meta$output.dir)
	return(thin.mcmc.set)
}

has.tfr.prediction <- function(mcmc=NULL, sim.dir=NULL) {
	if (!is.null(mcmc)) sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
	if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
	if(file.exists(file.path(sim.dir, 'predictions', 'prediction.rda'))) return(TRUE)
	return(FALSE)
}

get.tfr.prediction <- function(mcmc=NULL, sim.dir=NULL, mcmc.dir=NULL) {
	############
	# Returns an object of class bayesTFR.prediction
	# Set mcmc.dir to NA, if the prediction object should not have a pointer 
	# to the corresponding mcmc traces.
	############
	if (!is.null(mcmc)) 
		sim.dir <- if(is.character(mcmc)) mcmc else mcmc$meta$output.dir
	if (is.null(sim.dir)) stop('Either mcmc or directory must be given.')
	output.dir <- file.path(sim.dir, 'predictions')
	pred.file <- file.path(output.dir, 'prediction.rda')
	if(!file.exists(pred.file)) {
		warning('File ', pred.file, ' does not exist.')
		return(NULL)
	}
	load(file=pred.file)
	bayesTFR.prediction$output.directory <- output.dir
	
	pred <- bayesTFR.prediction
	# re-route mcmcs if necessary
	if(!is.null(mcmc.dir) || !has.tfr.mcmc(pred$mcmc.set$meta$output.dir)) {
		if((!is.null(mcmc.dir) && !is.na(mcmc.dir)) || is.null(mcmc.dir)) {
			new.path <- file.path(sim.dir, basename(pred$mcmc.set$meta$output.dir))
			if (has.tfr.mcmc(new.path)) pred$mcmc.set <- get.tfr.mcmc(new.path)
			else {
				est.dir <- if(is.null(mcmc.dir)) sim.dir else mcmc.dir
				pred$mcmc.set <- get.tfr.mcmc(est.dir)
			}
		}
	}
	return(pred)
}

.do.get.convergence.all <- function(type, package, sim.dir) {
	diag.dir <- file.path(sim.dir, 'diagnostics')
	if (!file.exists(diag.dir)) return(NULL)
	files <- list.files(diag.dir, pattern=paste('^', package, '[.]convergence_[0-9]+_[0-9]+[.]rda$', sep=''), 
							full.names=FALSE)
	if(length(files) <= 0) return(NULL)
	thin.matches <- regexpr('e_[0-9]+', files)
	bi.matches <- regexpr('[0-9]+[.]', files)
	result <- list()
	for(i in 1:length(files)) {
		thin <- as.integer(substr(files[i], start=thin.matches[i]+2, 
								stop=thin.matches[i]+attr(thin.matches, 'match.length')[i]-1))
		burnin <- as.integer(substr(files[i], start=bi.matches[i], 
								stop=bi.matches[i]+attr(bi.matches, 'match.length')[i]-2))
		result[[i]] <- do.call(paste('get.', type, '.convergence', sep=''), 
								list(sim.dir, thin=thin, burnin=burnin))
	}
	return(result)
}
get.tfr.convergence.all <- function(sim.dir=file.path(getwd(), 'bayesTFR.output')) {
	return(.do.get.convergence.all('tfr', 'bayesTFR', sim.dir=sim.dir))
}

get.tfr.convergence <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									thin=80, burnin=2000) {
	file.name <- file.path(sim.dir, 'diagnostics', paste('bayesTFR.convergence_', thin, '_', burnin, '.rda', sep=''))
	if(!file.exists(file.name)){
		warning('Convergence diagnostics in ', sim.dir, ' for burnin=', burnin, 
					' and thin=', thin, ' does not exist.')
		return(NULL)
	}
	bayesTFR.convergence <- local({
								load(file.name)
								bayesTFR.convergence})
	return(bayesTFR.convergence)
}

get.tfr3.convergence.all <- function(sim.dir=file.path(getwd(), 'bayesTFR.output')) {
	return(get.tfr.convergence.all(sim.dir=file.path(sim.dir, 'phaseIII')))
}

get.tfr3.convergence <- function(sim.dir=file.path(getwd(), 'bayesTFR.output'), 
									thin=60, burnin=10000) {
	return(get.tfr.convergence(file.path(sim.dir, 'phaseIII'), thin=thin, burnin=burnin))
}



has.mcmc.converged <- function(diag) return(diag$status['green'])	
get.burned.tfr.traces <- function(mcmc, par.names, burnin=0, thinning.index=NULL) {
	# get traces that are already loaded in the mcmc object
	traces <- mcmc$traces[, par.names]
	discard <- burnin - mcmc$traces.burnin
	if (discard > 0)
		traces <- traces[-seq(1, discard),]
	if(!is.null(thinning.index))
		traces <- traces[thinning.index,]
	return(traces)
}
"bdem.parameter.traces" <- function(mcmc, ...) UseMethod("bdem.parameter.traces")

bdem.parameter.traces.bayesTFR.mcmc <- function(mcmc, par.names, ...) {
	tran.names <- totran.names <- all.standard.names <- backtran.names <- tobacktran.names <- c()
	# Load traces from the disk
	if(is.null(mcmc$meta$phase) || mcmc$meta$phase == 2) {
		all.standard.names <- c(tfr.parameter.names(), get.trans.parameter.names(), 
							tfr.parameter.names.cs(), get.trans.parameter.names(cs=TRUE))
		tran.names <- c(get.trans.parameter.names(), get.trans.parameter.names(cs=TRUE))
		totran.names <- c(get.totrans.parameter.names(), get.totrans.parameter.names(cs=TRUE))
		backtran.names <- get.backtrans.parameter.names(cs=TRUE)
		tobacktran.names <- get.tobacktrans.parameter.names(cs=TRUE)
	} else {
		if(mcmc$meta$phase == 3) all.standard.names <- c(tfr3.parameter.names(), tfr3.parameter.names.cs())
	}
	return(.do.get.traces(mcmc, par.names=par.names, ..., 
							all.standard.names=all.standard.names,
							tran.names=tran.names, totran.names=totran.names,
							backtran.names=backtran.names, tobacktran.names=tobacktran.names))
}

.do.get.traces <- function(mcmc, par.names, file.postfix='', par.names.postfix='', burnin=0, 
							thinning.index=NULL, all.standard.names=c(), tran.names=c(), totran.names=c(), 
							backtran.names=c(), tobacktran.names=c()) {
	if (length(par.names) == 0) return (NULL)
	tran.names.l <- nchar(tran.names)
	ltran.names <- length(tran.names)
	totran.names.l <- nchar(totran.names)
	has.tran <- rep(FALSE, ltran.names)
	lbacktran.names <- length(backtran.names)
	backtran.names.l <- nchar(backtran.names)
	has.backtran <- rep(FALSE, lbacktran.names)
	values <- c()
 	valnames <- c()
 	loaded.files <- c()
 	compr.settings <- .get.compression.settings(mcmc$compression.type)
 	for(name in par.names) {
 		if(lbacktran.names > 0) {
 			new.has.backtran <- backtran.names %in% sapply(backtran.names.l, function(x) substr(name, 1, x))
 			has.backtran <- has.backtran | new.has.backtran
 			if(any(new.has.backtran)) next  # name is a back-transformation variable, therefore skip to next one
 		} 
 		if(ltran.names > 0) {
 			new.has.tran <- tran.names %in% sapply(tran.names.l, function(x) substr(name, 1, x))
 			has.tran <- has.tran | new.has.tran
 			if(any(new.has.tran)) next  # name is a transformation variable, therefore skip to next one
 		}
 		if (any(name == valnames)) next # name already loaded
 			
 		if (!any(name==all.standard.names)) {
 			#trim the name to the standard names, e.g. make 'gamma' from 'gamma_3'
 			name.to.load <- NULL
 			for (sname in all.standard.names) {
 				if (length(grep(paste('^', sname, '_', sep=''), name)) > 0) {
 					name.to.load <- sname
 					break
 				}
 			}
 			if (is.null(name.to.load)) {
 				warning('Parameter ', name, ' not found.')
 				next
 			}
 		} else {name.to.load <- name}
 		if (any(name.to.load == valnames)) next # name already loaded
 		#file.name <- paste(name.to.load, file.postfix, '.txt.bz2',sep='')
 		file.name <- paste(name.to.load, file.postfix, '.txt', compr.settings[2], sep='') 		
 		if (any(file.name == loaded.files)) next
 		#con <- bzfile(file.path(mcmc$meta$output.dir, mcmc$output.dir, file.name), open="rb")
 		#con <- file(file.path(mcmc$meta$output.dir, mcmc$output.dir, file.name), open="r")
 		con <- do.call(compr.settings[1], list(file.path(mcmc$meta$output.dir, mcmc$output.dir, file.name), 
 								open=paste("r", compr.settings[3], sep='')))
 		if(compr.settings[3]=='b')  # binary connection
 			raw.con <- textConnection(readLines(con))
 		else raw.con <- con
 		#close(con)
		vals <- as.matrix(read.table(file = raw.con))
		close(raw.con)
		#vals <- as.matrix(read.table(file = con))
		if(compr.settings[3]=='b') close(con)
		loaded.files <- c(loaded.files, file.name)
		if (dim(vals)[2] > 1) { #2d parameters get postfix
			valnames <- c(valnames, paste(name.to.load, '_', 1:dim(vals)[2], par.names.postfix, sep=''))
		} else {
			valnames <- c(valnames, paste(name.to.load, par.names.postfix, sep=''))
		}
		if (burnin > 0) {
			if (burnin > dim(vals)[1]) stop('Burnin is larger than the data size.')
			vals <- vals[-seq(1, burnin),,drop=FALSE]
		}
		if(!is.null(thinning.index))
			vals <- vals[thinning.index,,drop=FALSE]
		values <- cbind(values, vals)
	} # end of loading loop
	#stop('')
	colnames(values) <- valnames
	if(ltran.names == 0 && lbacktran.names == 0 ) return(values)
	# get transformed variables (alpha, delta, gamma)
	if(ltran.names > 0 && any(has.tran)) {
		for (i in which(has.tran)) {
			full.names <- grep(paste0(totran.names[i],'_'), valnames, value=TRUE)
			if (length(grep(paste0('^', totran.names[i], '(_|$)'), valnames)) == 0) { # load alpha/delta/gamma
				vals <- bdem.parameter.traces(mcmc, c(totran.names[i]), 
								file.postfix=file.postfix, par.names.postfix=par.names.postfix,
								burnin=burnin, thinning.index=thinning.index)
				full.names <- grep(paste0(totran.names[i],'_'), c(valnames, colnames(vals)), value=TRUE)
			} else {
				vals <- values[, full.names, drop=FALSE]
			}
			vals <- exp(vals)/apply(exp(vals),1,sum)
			values <- cbind(values, vals)
			valnames <- c(valnames, paste0(tran.names[i],'_', 
											substr(full.names, totran.names.l[i]+2, nchar(full.names))))
			colnames(values) <- valnames
		}
	}
	#stop('')
	# get back-transformed variables (Triangle_c1-3)
	if(lbacktran.names > 0 && any(has.backtran)) {
		vals <- c()
		for(i in 1:length(tobacktran.names)) {		
			full.names <- grep(paste0(tobacktran.names[i],'_'), valnames, value=TRUE)
			if (length(grep(paste0('^', tobacktran.names[i], '(_|$)'), valnames)) == 0) { # load gamma, U, Triangle_c4
				vals <- cbind(vals, bdem.parameter.traces(mcmc, c(tobacktran.names[i]), 
								file.postfix=file.postfix, par.names.postfix=par.names.postfix,
								burnin=burnin, thinning.index=thinning.index))
				full.names <- grep(paste0(tobacktran.names[i],'_'), c(valnames, colnames(vals)), value=TRUE)
			} else {
				vals <- cbind(vals, values[, full.names, drop=FALSE])
			}
		}
		vals <- .get.backtransformed.Triangles(vals, backtran.names[has.backtran])
		values <- cbind(values, vals)
	}
	return(values)
}

.get.backtransformed.Triangles <- function(data, par.names) {
	data.par.names <- colnames(data)
	# set names of the parameter columns for this country
	U.var <- grep('^U_c', data.par.names, value=TRUE) 
	d.var <- grep('^d_c', data.par.names, value=TRUE) 
	Triangle_c4.var <- grep("^Triangle_c4_c", data.par.names, value=TRUE) 
	gamma.vars <- sapply(paste0('^gammat_',1:3), function(x) grep(x, data.par.names, value=TRUE))
	# transform gamma_ci into Triangle_ci 
	# compute p_ci 
	#p_ci <- matrix(NA, nrow(data), 3) 
	#p_ci_denominator <- apply(exp(data[,gamma.vars]), 1, sum) 
	#for (i in 1:3) p_ci[,i] <- exp(data[,gamma.vars[i]])/p_ci_denominator 
	Delta <- data.frame(
    	Triangle_c1 = (data[,U.var] - data[,Triangle_c4.var])*data[,gamma.vars[1]],
    	Triangle_c2 = (data[,U.var] - data[,Triangle_c4.var])*data[,gamma.vars[2]],
    	Triangle_c3 = (data[,U.var] - data[,Triangle_c4.var])*data[,gamma.vars[3]]
    ) 
    Delta <- Delta[,par.names, drop=FALSE]
    suffix <- strsplit(data.par.names, "_")[[1]]
    suffix <- suffix[length(suffix)]
    colnames(Delta) <- paste(colnames(Delta), suffix, sep="_")
	return(as.matrix(Delta))
	
}


get.all.parameter.names <- function() {
    # First element in each tuple indicates if the parameter is transformable, 
    # the second how many values it has.
	return(list(alpha=c(TRUE, 3),  delta=c(FALSE, 3), Triangle4=c(FALSE, 1), delta4=c(FALSE, 1), 
				psi=c(FALSE, 1), chi=c(FALSE, 1), a_sd=c(FALSE, 1), b_sd=c(FALSE, 1), 
				const_sd=c(FALSE, 1), S_sd=c(FALSE, 1), sigma0=c(FALSE, 1), mean_eps_tau=c(FALSE, 1),
				sd_eps_tau=c(FALSE, 1)))
}				

get.all.parameter.names.cs <- function() {
    # First element in each tuple indicates if the parameter is transformable, 
    # the second how many values it has.
	return(list(gamma=c(1, 3), U=c(0, 1), d=c(0, 1), Triangle_c4=c(0, 1), 
				Triangle_c1=c(2, 1), Triangle_c2=c(2, 1), Triangle_c3=c(2, 1)))
}

get.totrans.parameter.names <- function(cs=FALSE) {
	pars <- if(cs) get.all.parameter.names.cs() else get.all.parameter.names() 
	is.trans <- sapply(pars, function(x) return(x[1] == 1))
	return(names(pars)[is.trans])
}

get.trans.parameter.names <- function(cs=FALSE) return(paste(get.totrans.parameter.names(cs), 't', sep=''))

get.backtrans.parameter.names <- function(cs=FALSE) {
	if(!cs) return (c())
	pars <- get.all.parameter.names.cs()
	is.backtrans <- sapply(pars, function(x) return(x[1] == 2))
	return(names(pars)[is.backtrans])
}

get.tobacktrans.parameter.names <- function(cs=FALSE) {
	if(!cs) return (c())
	return(c("U", "gammat", "Triangle_c4"))
}

get.other.parameter.names <- function(cs=FALSE) {
	pars <- if(cs) get.all.parameter.names.cs() else get.all.parameter.names()
	is.trans <- sapply(pars, function(x) return(x[1]) > 0)
	return(names(pars)[!is.trans])
}

get.all.parameter.names.extended <- function(cs=FALSE) {
	pars <- c()
	all.pars <- if(cs) get.all.parameter.names.cs() else get.all.parameter.names()
	for (ipar in 1:length(all.pars)) {
		par <- all.pars[ipar]
		paropt <- unlist(par)
		name <- names(par)
		pars <- c(pars, if (paropt[2] > 1) paste(name, 1:paropt[2], sep='_') else name)
		if (paropt[1]) {
			name <- paste(name,'t', sep='')
			pars <- c(pars, if (paropt[2] > 1) paste(name, 1:paropt[2], sep='_') else name)
		}
	}
	return(pars)
}

tfr.parameter.names.extended <- function() {
	# return a list of all parameters with those parameters extended by i
	# that have more than 1 value, e.g. alpha_2, alphat_1, delta_3
	return(get.all.parameter.names.extended(cs=FALSE))
}

tfr.parameter.names.cs.extended <- function(country.code=NULL) {
	# return a list of cs-parameters with gamma extended by i
	pars <- get.all.parameter.names.extended(cs=TRUE)
	if(!is.null(country.code)) pars <- paste(pars, '_c', country.code, sep='')
	return(pars)
}

tfr.parameter.names <- function(trans=NULL) {
	# Return all country-independent parameter names. 
	# trans can be NULL or logical.
	# If 'trans' is TRUE,
	# names of the transformable parameters (alpha) are replaced by the names 
	# of the transformed parameters.
	# If 'trans' is NULL, all parameter names, 
	# including the transformable parameters are returned.
	other.pars <- get.other.parameter.names(cs=FALSE)
	if (is.null(trans)) return (c(get.totrans.parameter.names(), get.trans.parameter.names(), other.pars))
	if (trans) return (c(get.trans.parameter.names(), other.pars))
	return (c(get.totrans.parameter.names(), other.pars))
}
	
tfr.parameter.names.cs <- function(country.code=NULL, trans=NULL, back.trans=TRUE) {
	#Return all country-specific parameter names. 
	# See comments in tfr.parameter.names(). Transformable parameter is gamma.
	#If country is not NULL, it must be a country code.
	#It is attached to the parameter name.
	par.names <- get.other.parameter.names(cs=TRUE)
	backtrans.names <- if(back.trans) get.backtrans.parameter.names(cs=TRUE) else c()
	if (is.null(trans)) par.names <- c(par.names, backtrans.names,
										get.totrans.parameter.names(cs=TRUE), 										
										get.trans.parameter.names(cs=TRUE))
	else par.names <- c(par.names, backtrans.names,
								if (trans) get.trans.parameter.names(cs=TRUE) 
									else get.totrans.parameter.names(cs=TRUE))
	if (is.null(country.code))
		return(par.names)
	return(paste(par.names, '_c', country.code, sep=''))
}

tfr3.parameter.names <- function() return(c('mu', 'rho', 'sigma.mu', 'sigma.rho', 'sigma.eps'))
tfr3.parameter.names.cs <- function(country.code=NULL) return(paste(c('mu', 'rho'), '.c', country.code,sep=''))

get.total.iterations <- function(mcmc.list, burnin=0) {
	# Return total number of iterations sum-up over chains after discarding burnin in each chain
	get.iter <- function(x) return(x$finished.iter - burnin)
	return(sum(sapply(mcmc.list, get.iter)))
}

get.thinned.burnin <- function(mcmc, burnin) {
	if (burnin==0) return(0)
	if (mcmc$thin == 1) return(burnin)
	return(1 + if(burnin >= mcmc$thin) length(seq(mcmc$thin, burnin, by=mcmc$thin)) else 0)
}

get.stored.mcmc.length <- function(mcmc.list, burnin=0) {
	# Return total number of iterations sum-up over chains after discarding burnin in each chain,
	# taking into account the original value of thin.
	# It should correspond to the total length of all chains stored on disk, minus burnin.
	get.iter <- function(x) return(x$length - get.thinned.burnin(x, burnin))
	return(sum(sapply(mcmc.list, get.iter)))
}

do.get.tfr.parameter.traces <- function(is.cs, mcmc.list, par.names, country.obj=NULL, 
										burnin=0, thinning.index=NULL, thin=NULL) {
	# get parameter traces either from disk or from memory (if they were already loaded)
	# par.names are either country-independent (if is.cs is FALSE), or country-specific (if is.cs is TRUE)
	values <- c()
	if (is.null(thinning.index) && !is.null(thin) && thin > mcmc.list[[1]]$thin) {
		total.iter <- get.stored.mcmc.length(mcmc.list, burnin)
		thinning.index <- unique(round(seq(1, total.iter, by=thin/mcmc.list[[1]]$thin)))
	}
	at.iter <- 1
	for(mcmc in mcmc.list) {
		this.thinning.index <- NULL
		th.burnin <- get.thinned.burnin(mcmc, burnin)
		if(!is.null(thinning.index)) {
			this.thinning.index <- thinning.index[(thinning.index >= at.iter) & 
									(thinning.index < at.iter+mcmc$length-th.burnin)] - at.iter+1
			if (length(this.thinning.index) == 0) {
				at.iter <- at.iter+mcmc$length-th.burnin
				next
			}
		}
    	if (no.traces.loaded(mcmc) || th.burnin < mcmc$traces.burnin) {
    		traces <- if(is.cs) load.tfr.parameter.traces.cs(mcmc, country.obj$code, par.names, burnin=th.burnin, 
    											thinning.index=this.thinning.index)
    					else bdem.parameter.traces(mcmc, par.names, burnin=th.burnin, thinning.index=this.thinning.index)
        } else {
            traces <- if(is.cs) get.burned.tfr.traces(mcmc, get.full.par.names.cs(par.names, colnames(mcmc$traces), 
            									country=country.obj$index), 
            								th.burnin, thinning.index=this.thinning.index)
            		  else get.burned.tfr.traces(mcmc, par.names, th.burnin, thinning.index=this.thinning.index)
        }
       	values <- rbind(values, traces)
       	at.iter <- at.iter+mcmc$length-th.burnin
    }
    return(values)
}

get.tfr.parameter.traces <- function(mcmc.list, par.names=tfr.parameter.names(), 
										burnin=0, thinning.index=NULL, thin=NULL) {
	# get parameter traces either from disk or from memory, if they were already loaded
	return(do.get.tfr.parameter.traces(is.cs=FALSE, mcmc.list=mcmc.list, par.names=par.names, 
										burnin=burnin, thinning.index=thinning.index, thin=thin))
}

get.tfr.parameter.traces.cs <- function(mcmc.list, country.obj, par.names=tfr.parameter.names.cs(), 
										burnin=0, thinning.index=NULL, thin=NULL) {
	# country.obj is result of get.country.object()
	# get traces for country-specific parameters either from disk or from memory, if they were already loaded
	return(do.get.tfr.parameter.traces(is.cs=TRUE, mcmc.list=mcmc.list, par.names=par.names, country.obj=country.obj,
										burnin=burnin, thinning.index=thinning.index, thin=thin))
}

get.tfr3.parameter.traces <- function(mcmc.list, par.names=tfr3.parameter.names(), ...)
	return(get.tfr.parameter.traces(mcmc.list, par.names, ...))
	
get.tfr3.parameter.traces.cs <- function(mcmc.list, country.obj, par.names=tfr3.parameter.names.cs(), ...)
	return(get.tfr.parameter.traces.cs(mcmc.list, country.obj, par.names, ...))

load.tfr.parameter.traces <- function(mcmc, par.names=tfr.parameter.names(), burnin=0, thinning.index=NULL) 
 	return(bdem.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index))

load.tfr.parameter.traces.cs <- function(mcmc, country, par.names=tfr.parameter.names.cs(), burnin=0, 
										thinning.index=NULL) {
 	return(bdem.parameter.traces(mcmc, par.names, paste('_country', country, sep=''),
						par.names.postfix=paste('_c', country, sep=''), burnin=burnin, 
						thinning.index=thinning.index))
}

load.tfr3.parameter.traces <- function(mcmc, par.names=tfr3.parameter.names(), ...)
	return(load.tfr.parameter.traces(mcmc, par.names=par.names, ...))

load.tfr3.parameter.traces.cs <- function(mcmc, country, par.names=tfr3.parameter.names.cs(), ...)
	return(load.tfr.parameter.traces.cs(mcmc, country=country, par.names=par.names, ...))

load.tfr.parameter.traces.all <- function(mcmc, par.names=tfr.parameter.names(), 
										 par.names.cs=tfr.parameter.names.cs(),
										 burnin=0, thinning.index=NULL) {
	result <- load.tfr.parameter.traces(mcmc, par.names, burnin=burnin, thinning.index=thinning.index)
	if(!is.null(par.names.cs))
		for (country in get.countries.index(mcmc$meta)) {
			result <- cbind(result, 
						load.tfr.parameter.traces.cs(mcmc, 
												    get.country.object(country, 
												         mcmc$meta, index=TRUE)$code, 
												    par.names.cs, burnin=burnin,
												    thinning.index=thinning.index))
		}
	return (result)
}

load.tfr3.parameter.traces.all <- function(mcmc, par.names=tfr3.parameter.names(), 
										 par.names.cs=tfr3.parameter.names.cs(), ...)
	return(load.tfr.parameter.traces.all(mcmc, par.names=par.names, par.names.cs=par.names.cs, ...))
	
get.full.par.names.cs <- function(par.names, full.par.names, country=NULL, index=FALSE) {
	# Return full name of par.names that are included in full.par.names
	# which are suppose to be all country-specific parameters.
	# E.g. for 'gamma', it would return 'gamma_1_c1', ..., gamma_1_cX'
	# If index is TRUE, return the index of the matches.
	result <- c()	
	for (name in par.names) {
		pattern <- paste('^', name, '_.*c', sep='')
		if (!is.null(country)) {
			pattern <- paste(pattern, country, '$', sep='')
		}
		result <- c(result, grep(pattern, full.par.names, value=!index))
	}
	return(result)
}

get.full.par.names <- function(par.names, full.par.names, index=FALSE) {
	# Return full name of par.names that are included in full.par.names
	# which are suppose to be all country-independent parameters.
	# E.g. for 'alpha', it would return 'alpha_1', alpha_2, alpha_3'
	# If index is TRUE, return the index of the matches.
	result <- c()	
	for (name in par.names) {
		pattern <- paste('^', name, '(_|$)', sep='') # either has '_' suffix or matches exactly
		result <- c(result, grep(pattern, full.par.names, value=!index))
	}
	return(result)
}

"coda.mcmc" <- function(mcmc, ...) UseMethod("coda.mcmc")

coda.mcmc.bayesTFR.mcmc <- function(mcmc, country=NULL, par.names=NULL, 
						par.names.cs=NULL, burnin=0, thin=1, ...
						) {
	# Return a coda object for this mcmc and parameter names
	index <- NULL
	btobject <- burn.and.thin(mcmc, burnin=burnin, thin=thin)
	thin <- btobject$thin
	th.burnin <- btobject$burnin
	if(!is.null(btobject$index)) index <- btobject$index
	if(missing(par.names)) 
		par.names <- if(!is.null(mcmc$meta$phase) && mcmc$meta$phase == 3)
			tfr3.parameter.names() else tfr.parameter.names()
	if(missing(par.names.cs)) 
		par.names.cs <- if(!is.null(mcmc$meta$phase) && mcmc$meta$phase == 3)
			tfr3.parameter.names.cs() else tfr.parameter.names.cs()

	if (!is.null(country)) { # for specific country
		if (burnin < mcmc$traces.burnin || no.traces.loaded(mcmc)) {
			values <- load.tfr.parameter.traces(mcmc, par.names, burnin=th.burnin, thinning.index=index)
			values <- cbind(values, load.tfr.parameter.traces.cs(mcmc, 
								get.country.object(country, mcmc$meta)$code, par.names.cs, 
								burnin=th.burnin, thinning.index=index))
		} else {
			values <- get.burned.tfr.traces(mcmc, 
							c(get.full.par.names(par.names, colnames(mcmc$traces)), 
							  get.full.par.names.cs(par.names.cs, colnames(mcmc$traces), 
							  	get.country.object(country, mcmc$meta)$code)), th.burnin,
							  	thinning.index=index)
		}
	} else { #no country specified
		if (no.traces.loaded(mcmc)) { # traces not loaded
			values <- load.tfr.parameter.traces.all(mcmc, par.names, par.names.cs, burnin=th.burnin, 
												    thinning.index=index)
		} else { # traces loaded but get the right burnin
			values <- get.burned.tfr.traces(mcmc, c(get.full.par.names(par.names, colnames(mcmc$traces)), 
                             get.full.par.names.cs(par.names.cs, colnames(mcmc$traces))), th.burnin,
                             thinning.index=index)
        }
    }
    # filter out unnecessary parameters
    values <- filter.traces(values, c(par.names, par.names.cs))
    return(mcmc(values, end=mcmc$finished.iter, thin=thin, ...))
}
	
	
coda.list.mcmc <- function(mcmc=NULL, country=NULL, chain.ids=NULL,
							sim.dir=file.path(getwd(), 'bayesTFR.output'), 
							par.names=tfr.parameter.names(), 
							par.names.cs=tfr.parameter.names.cs(), 
							rm.const.pars=FALSE,
							burnin=0, 
							low.memory=FALSE, ...) {
	# return a list of mcmc objects that can be analyzed using the coda package
	mcmc.list <- mcmc
	if (is.null(mcmc.list)) {
		mcmc.list <- get.tfr.mcmc(sim.dir, chain.ids=chain.ids, low.memory=low.memory, 
									burnin=burnin)$mcmc.list
	} else {
		mcmc.list <- get.mcmc.list(mcmc.list)
		if (!is.null(chain.ids)) {
			mcmc.list <- mcmc.list[chain.ids]
		}
	}
	result <- list()
	i <- 1
	for(mc in mcmc.list) {
		result[[i]] <- coda.mcmc(mc, country=country, par.names=par.names, par.names.cs=par.names.cs,
								 burnin=burnin, ...)

		if (rm.const.pars) {
			# remove parameters that are constant for all iterations
			if (is.null(dim(result[[i]]))) { # one parameter in the chain
				if(all(result[[i]]==result[[i]][1]))  { # no parameters are kept
					result[[i]] <- NULL
					i <- i - 1
				}
			} else { # multiple parameters in the chain
				if (dim(result[[i]])[2]==1) {
					colindex <- !all(result[[i]][,1] == result[[i]][1,1])
				} else {
					colindex <- !apply(t(apply(result[[i]], 1, '==', result[[i]][1,])), 2, all)
				}
				if (sum(colindex) == 0) { # no parameters are kept
					result[[i]] <- NULL
					i <- i - 1
				} else {
					result[[i]] <- result[[i]][,colindex]
				}
			}
		}
		i <- i+1
	}
	return(mcmc.list(result))
}

coda.list.mcmc3 <- function(mcmc=NULL, country=NULL, chain.ids=NULL,
							sim.dir=file.path(getwd(), 'bayesTFR.output'), 
							par.names=tfr3.parameter.names(), 
							par.names.cs=tfr3.parameter.names.cs(), 
							burnin=0, low.memory=FALSE, ...) {
	if (is.null(mcmc)) 
		mcmc <- get.tfr3.mcmc(sim.dir, chain.ids=chain.ids, low.memory=low.memory, 
									burnin=burnin)$mcmc.list
	else
		if(class(mcmc)=='bayesTFR.prediction')
			stop('Function not available for bayesTFR.prediction objects.')
	return(coda.list.mcmc(mcmc=mcmc, country=country, chain.ids=chain.ids, sim.dir=NULL, 
			par.names=par.names, par.names.cs=par.names.cs, rm.const.pars=FALSE, burnin=burnin, ...))										
}

filter.traces <- function(values, par.names) {
	valuenames <- colnames(values)
    lpar.names <- length(par.names)
    include <- rep(FALSE, length(valuenames))
    for(i in 1:lpar.names) {
    	idx <- grep(paste('^', par.names[i], '$', sep=''), valuenames) # exact match
    	if (length(idx) > 0) {
    		include[idx] <- TRUE
    		next
    	}
    	idx <- grep(paste0('^', par.names[i], '_'), valuenames) # partial match
    	if (length(idx) > 0) include[idx] <- TRUE
    }
	#v <- array(values[,include], c(nrow(values), sum(include)))
	v <- values[,which(include), drop=FALSE]
	#colnames(v) <- colnames(values)[include]
    return(v)
}

"get.mcmc.list" <- function(mcmc.list, ...) UseMethod("get.mcmc.list")

get.mcmc.list.bayesTFR.mcmc.set <- function(mcmc.list, ...) return(mcmc.list$mcmc.list)
get.mcmc.list.bayesTFR.mcmc <- function(mcmc.list, ...) return(list(mcmc.list))
get.mcmc.list.bayesTFR.prediction <- function(mcmc.list, ...) return(mcmc.list$mcmc.set$mcmc.list)
get.mcmc.list.list <- function(mcmc.list, ...) return(mcmc.list)

"get.mcmc.meta" <- function(meta, ...) UseMethod("get.mcmc.meta")
get.mcmc.meta.bayesTFR.mcmc.set <- function(meta, ...) return(meta$meta)
get.mcmc.meta.bayesTFR.mcmc.meta <- function(meta, ...) return(meta)
get.mcmc.meta.bayesTFR.mcmc <- function(meta, ...) return(meta$meta)
get.mcmc.meta.bayesTFR.prediction <- function(meta, ...) return(meta$mcmc.set$meta)


get.country.object <- function(country, meta=NULL, country.table=NULL, index=FALSE) {
	# If meta object is not given, country.table containing columns 'code' and 'name' must be given. 
	# If 'country' is numeric, 'index' determines if 'country' is an index (TRUE) or code (FALSE)
	if (!is.null(meta)) {
		codes <- meta$regions$country_code
		names <- meta$regions$country_name
	} else {
		codes <- country.table[,'code']
		names <- country.table[,'name']	
	}
	l <- length(codes)
	found <- TRUE
	if (is.numeric(country)) {
		if (index) {
			country.idx <- country
			country.code <- codes[country.idx]
		} else { 
			country.code <- country
			country.idx <- (1:l)[codes==country.code]
		}
		country.name <- as.character(names[country.idx])
		if (length(country.name) == 0) found <- FALSE
	} else {
		country.name <- country
		country.idx <- (1:l)[names==country.name]
		country.code <- codes[country.idx]
		if (length(country.idx) == 0) found <- FALSE
	}
	if (!found) 
		country.name <- country.idx <- country.code <- NULL
	return(list(name=country.name, index=country.idx, code=country.code))
}

country.names <- function(meta, countries=NULL, index=FALSE) {
	meta <- get.mcmc.meta(meta)
	if (is.null(countries)) return(meta$regions$country_name)
	if (index) return(meta$regions$country_name[countries])
	get.country.index <- function(code) {
		return(get.country.object(code, meta)$index)
		}
	return(meta$regions$country_name[sapply(countries, get.country.index)])
}

summary.bayesTFR.mcmc <- function(object, country=NULL, 
								par.names=NULL, par.names.cs=NULL, 
								thin=1, burnin=0, ...) {
	if(is.null(country) && missing(par.names.cs)) par.names.cs <- NULL
	if(is.null(par.names))
	 	par.names <- if(is.null(object$meta$phase) || object$meta$phase == 2) tfr.parameter.names(trans=TRUE) 
	 					else tfr3.parameter.names()
	if(is.null(par.names.cs) && !is.null(country))
		par.names.cs <- if(is.null(object$meta$phase) || object$meta$phase == 2) tfr.parameter.names.cs(trans=TRUE) 
	 					else tfr3.parameter.names.cs()
	if (!is.null(country)) {
		country.obj <- get.country.object(country, object$meta)
		cat('\nCountry:', country.obj$name, '\n')
		if (!is.element(country.obj$index, object$meta$id_DL)) {
			cat('\tnot used for estimation because no decline observed.\n')
			return(NULL)
		}
		country <- country.obj$code
	} 
	summary(coda.mcmc(object, country=country, par.names=par.names,
							par.names.cs=par.names.cs, thin=thin, burnin=burnin), ...)
}

summary.bayesTFR.mcmc.set <- function(object, country=NULL, chain.id=NULL, 
								par.names=NULL, par.names.cs=NULL, 
								meta.only=FALSE, thin=1, burnin=0, ...) {
	if(is.null(country) && missing(par.names.cs)) par.names.cs <- NULL
	if(is.null(object$meta$phase) || object$meta$phase == 2) {
		if(is.null(par.names)) par.names <- tfr.parameter.names(trans=TRUE)
		if(!is.null(country) && is.null(par.names.cs)) par.names.cs <- tfr.parameter.names.cs(trans=TRUE)
		.summary.mcmc.set.phaseII(object, country, chain.id, par.names, par.names.cs, meta.only, thin, burnin, ...)
	} else { # phase III
		if(is.null(par.names)) par.names <- tfr3.parameter.names()
		if(!is.null(country) && is.null(par.names.cs)) par.names.cs <- tfr3.parameter.names.cs()
		.summary.mcmc.set.phaseIII(object, country, chain.id, par.names, par.names.cs, meta.only, thin, burnin, ...)
	}
}

.summary.mcmc.set.phaseII <- function(object, country=NULL, chain.id=NULL, 
								par.names=tfr.parameter.names(trans=TRUE), 
								par.names.cs=tfr.parameter.names.cs(trans=TRUE), 
								meta.only=FALSE, thin=1, burnin=0, ...) {
	cat('\nMCMCs of phase II')
	cat('\n=================')
	cat('\nNumber of countries:', object$meta$nr_countries)
	cat('\nHyperparameters estimated using', 
		length(object$meta$id_DL[object$meta$id_DL<=object$meta$nr_countries_estimation]), 
			'countries.\n')
	cat('\nWPP:', object$meta$wpp.year)
	cat('\nInput data: TFR for period', object$meta$start.year, '-', object$meta$present.year,'.')
	cat('\n')
	if(meta.only) {
		return(get.meta.only(object))
	} 
	if (!is.null(chain.id))
		return(summary(object$mcmc.list[[chain.id]], country=country, par.names=par.names,
							par.names.cs=par.names.cs, thin=thin, burnin=burnin, ...))
	if (!is.null(country)) {
		country.obj <- get.country.object(country, object$meta)
		cat('\nCountry:', country.obj$name, '\n')
		if (!is.element(country.obj$index, object$meta$id_DL)) {
			cat('\tnot used for estimation because no decline observed.\n')
			return(NULL)
		}
		country <- country.obj$code
	}
	summary(coda.list.mcmc(object, country=country, par.names=par.names,
							par.names.cs=par.names.cs, thin=thin, burnin=burnin), ...)		
}

.summary.mcmc.set.phaseIII <- function(object, country=NULL, chain.id=NULL, 
								par.names=NULL, par.names.cs=NULL, 
								meta.only=FALSE, thin=1, burnin=0, ...) {
	cat('\nMCMCs of phase III')
	cat('\n==================')
	cat('\nNumber of countries:', object$meta$nr.countries)
	cat('\nNumber of observations:', sum(sapply(object$mcmc.list[[1]]$observations, function(x) length(x)-1)))
	cat('\nWPP:', object$meta$parent$wpp.year)
	cat('\n')
	if(meta.only) return(get.meta.only(object))
	if (!is.null(chain.id))
		return(summary(object$mcmc.list[[chain.id]], country=country, par.names=par.names,
							par.names.cs=par.names.cs, thin=thin, burnin=burnin, ...))
	if (!is.null(country)) {
		country.obj <- get.country.object(country, object$meta)
		cat('\nCountry:', country.obj$name, '\n')
		if (!is.element(country.obj$index, object$meta$id_phase3)) {
			cat('\tnot used for estimation because it has not reached phase III yet.\n')
			return(NULL)
		}
		country <- country.obj$code
	}
	summary(coda.list.mcmc(object, country=country, par.names=par.names,
							par.names.cs=par.names.cs, thin=thin, burnin=burnin), ...)	
}


get.meta.only <- function(object) {
	get.iter <- function(x) x$finished.iter
	res <- list(nr.chains=object$meta$nr.chains, 
					iters=sapply(object$mcmc.list, get.iter),
					thin=object$mcmc.list[[1]]$thin)
	class(res) <- paste('summary.', class(object), '.meta', sep='')
	return(res)
}
 
print.summary.bayesTFR.mcmc.set.meta <- function(x, ...) {
	cat('\nNumber of chains =', x$nr.chains)
	cat('\nIterations =', paste(1,':',sum(x$iters)))
	cat('\nThinning interval =', x$thin)
	cat('\nChains sample sizes:', paste(x$iters, collapse=', '))
	cat('\n')
	invisible(x)
}

get.prediction.summary.data <- function(object, unchanged.pars, country, compact) {
	res <- list()
	for (par in unchanged.pars)
		res[[par]] <- object[[par]]
	proj.and.present.years <- if(is.null(object$proj.years)) 
				seq(object$end.year-5*object$nr.projections-2, 
										object$end.year-2, by=5)
							else object$proj.years

	res$projection.years <- proj.and.present.years[2:length(proj.and.present.years)]
	if(is.null(country)) return(res)
	if (!is.list(country))
		country <- get.country.object(country, object$mcmc.set$meta)
	if(is.null(country$code)) stop('No prediction available for this country.')
	res$country.name <- country$name
	
	if(compact) {
		quantiles.to.show <- c(0.025, 0.05, 0.1, 0.25, 0.5, 0.75, 0.9, 0.95, 0.975)
		quant.index <- is.element(dimnames(object$quantiles)[[2]], as.character(quantiles.to.show))
	} else {
		quant.index <- 1:dim(object$quantiles)[2]
	}
	res$projections <- cbind(t(object$traj.mean.sd[country$index,,]), t(object$quantiles[country$index,quant.index,]))
	colnames(res$projections) <- c('mean', 'SD', 
			paste(as.numeric(dimnames(object$quantiles)[[2]][quant.index])*100, '%', sep=''))
	rownames(res$projections) <- proj.and.present.years
	return(res)
}

.update.summary.data.by.shift <- function(res, object, country) {
	if(is.null(res$projections)) return(res)
	country <- get.country.object(country, object$mcmc.set$meta)
	shift <- get.tfr.shift(country$code, object)
	res$manual <- FALSE
	if(!is.null(shift)) {
		res$projections[,c(1,3:dim(res$projections)[2])] <- res$projections[,c(1,3:dim(res$projections)[2])] + t(matrix(shift,
												 ncol=nrow(res$projections), nrow=ncol(res$projections)-1, byrow=TRUE))
		res$manual <- TRUE
	}
	return(res)
}

summary.bayesTFR.prediction <- function(object, country=NULL, compact=TRUE, ...) {
	res <- get.prediction.summary.data(object, 
				unchanged.pars=c('burnin', 'thin', 'nr.traj', 'mu', 'rho', 'sigmaAR1', 'use.tfr3', 'burnin3', 'thin3'), 
				country=country, compact=compact)
	class(res) <- 'summary.bayesTFR.prediction'
	return(.update.summary.data.by.shift(res, object, country))
}

print.summary.bayesTFR.prediction <- function(x, digits = 3, ...) {
	cat('\nProjections:', length(x$projection.years), '(', x$projection.years[1], '-', 
					x$projection.years[length(x$projection.years)], ')')
	cat('\nTrajectories:', x$nr.traj)
	cat('\nPhase II burnin:', x$burnin)
	cat('\nPhase II thin:', x$thin)
	if(x$use.tfr3) {
		cat('\nPhase III burnin:', x$burnin3)
		cat('\nPhase III thin:', x$thin3)
	} else {
		cat('\nParameters of AR(1):\n')
		arpars <- c(x$mu, x$rho, x$sigmaAR1)
		ar.table<-matrix(arpars, nrow=1, ncol=length(arpars),
						dimnames=list('',c('mu', 'rho', rep('sigma', length(x$sigmaAR1)))))
		print(ar.table, digits=digits, ...)
	}
	if(!is.null(x$country.name)) {
		cat('\nCountry:', x$country.name, '\n')
		cat('\nProjected TFR:')
		if(x$manual) cat(' (values have been manually modified)')
		cat('\n')
		print(x$projections, digits=digits, ...)
	}
}


summary.bayesTFR.convergence <- function(object, expand=FALSE, ...) {
	cat('\nConvergence diagnostics for burnin =', object$burnin, 'and thin =', object$thin)
	cat('\n********************************************************')
	cat('\nFull chains:')
	cat('\n============')
	print(summary(object$mcmc.set, meta.only=TRUE))

	if(!is.null(object$thin.mcmc)) {
		cat('\nThinned, burned and collapsed chain:')
		cat('\n====================================')
		print(get.meta.only(object$thin.mcmc))
	}
	
	if(object$express) cat('\nExpress diagnostics - no country-specific parameters were included.')
	else {
		cat('\nConvergence checked on', object$nr.countries['used'], 
					'countries out of', object$nr.countries['total'], 'countries total.')
	}
	cat('\n')
	if(nrow(object$result) > 0) {
		if (object$lresult.country.independent > 0) {
			cat('\nNot converged country-independent parameters:')
			cat('\n---------------------------------------------\n')
			print(object$result[1:object$lresult.country.independent,,drop=FALSE])
		}
		if (nrow(object$result) - object$lresult.country.independent > 0) {
			cat('\nNot converged country-specific parameters:')
			cat('\n-------------------------------------------\n')
			print(object$result[(object$lresult.country.independent+1):nrow(object$result),,drop=FALSE])
		}
		cat('\n\nAt least', object$iter.needed, 'more iterations needed  to achieve convergence.\n')
	}
	if(object$status['green']) {
		cat('\nSimulation has converged.')
		cat('\nNumber of trajectories to be used:', object$use.nr.traj)
	}
	if(expand && !object$status['green'] && !is.null(object$country.specific) 
			&& !is.null(object$country.specific$not.converged.parameters)) {
		cat('\nWarning: ')
		cat('The following parameters did not converge:\n')
		print(object$country.specific$not.converged.parameters)
	}		
	cat('\nStatus:', names(object$status)[object$status], '\n')
}

tfr.info <- function(sim.dir) {
	mc <- get.tfr.mcmc(sim.dir=sim.dir)
	if (is.null(mc)) {
		cat('No TFR simulation available in', sim.dir)
		return (NULL)
	}
	summary(mc)
}

get.cov.gammas <- function(mcmc.set=NULL, sim.dir = NULL, burnin = 200, chain.id=1){
	# this is for one chain
	if (is.null(mcmc.set))
		mcmc.set <- get.tfr.mcmc(sim.dir=sim.dir)
 	cov_gammas_cii = array(NA, c(mcmc.set$meta$nr_countries, 3,3))
 	for (country in mcmc.set$meta$id_DL){
 		country.obj <- get.country.object(country, mcmc.set$meta, index=TRUE)
  		gamma_si <- load.tfr.parameter.traces.cs(mcmc.set$mcmc.list[[chain.id]], 
  							par.names='gamma', country=country.obj$code, burnin=burnin)
        cov_gammas_cii[country,, ] <- cov(gamma_si)
 	}
 	cov_gammas_cii = (2.4^2/3)*cov_gammas_cii 
 	return(list(values=cov_gammas_cii, country_codes=mcmc.set$meta$regions$country_code))
}

get.thinning.index <- function(nr.points, all.points) {
 	if (!is.null(nr.points)) {
		nr.points <- ifelse(nr.points >= all.points, all.points, nr.points)
	} else {
		nr.points <- all.points
	}
	if (nr.points > 0) {
		step <- all.points/nr.points
		idx <- floor(seq(floor(step), all.points, by=step))
	} else idx<-NULL
	return(list(nr.points=nr.points, index=idx))
}

burn.and.thin <- function(mcmc, burnin=0, thin=1) {
	# Return thin and burnin that is consolidated with the original thin usd for storing mcmc.
	# If there is need for more thinning, it returns the corresponding thinning index
	th.burnin <- get.thinned.burnin(mcmc,burnin)
	index <- NULL
	if (thin > mcmc$thin) {
    	index <- unique(round(seq(1,mcmc$length-th.burnin, by=thin/mcmc$thin)))
    }
    thin <- max(mcmc$thin, thin) # thin cannot be smaller than the original thin
	return(list(thin=thin, burnin=th.burnin, index=index))
}

no.traces.loaded <- function(mcmc) return((length(mcmc$traces) == 1) && mcmc$traces == 0)

tfr.set.identical <- function(mcmc.set1, mcmc.set2, include.output.dir=TRUE) {
	# Test if two bayesTFR sets are identical
	if(!include.output.dir) 
		mcmc.set1$meta$output.dir <- mcmc.set2$meta$output.dir <- NULL
	same <- setequal(names(mcmc.set1), names(mcmc.set2)) && identical(mcmc.set1$meta, mcmc.set2$meta) && length(mcmc.set1$mcmc.list) == length(mcmc.set2$mcmc.list)
	if(!same) return(same)
	for(i in 1:length(mcmc.set1$mcmc.list)) {
		if(!include.output.dir) mcmc.set1$mcmc.list[[i]]$meta$output.dir <- mcmc.set2$mcmc.list[[i]]$meta$output.dir <- NULL
		same <- same && tfr.identical(mcmc.set1$mcmc.list[[i]], mcmc.set2$mcmc.list[[i]])
	}
	return(same)
}

tfr.identical <- function(mcmc1, mcmc2) {
	# Test if two mcmcs are identical	
	same <- setequal(names(mcmc1), names(mcmc2))
	for(item in names(mcmc1)) 
		same <- same && identical(mcmc1[[item]], mcmc2[[item]])
	return(same)
}

get.tfr.trajectories <- function(tfr.pred, country) {
	country.obj <- get.country.object(country, tfr.pred$mcmc.set$meta)
	return(get.trajectories(tfr.pred, country.obj$code)$trajectories)
}

"get.nr.countries" <- function(meta, ...) UseMethod("get.nr.countries")
 
get.nr.countries.bayesTFR.mcmc.meta <- function(meta, ...) 
	return (if(is.null(meta$phase) || (meta$phase==2)) meta$nr_countries else meta$nr.countries)

"get.nr.countries.est" <- function(meta, ...) UseMethod("get.nr.countries.est")
 
get.nr.countries.est.bayesTFR.mcmc.meta <- function(meta, ...) 
	return (if(is.null(meta$phase) || (meta$phase==2)) meta$nr_countries_estimation else meta$nr.countries)

"get.data.matrix" <- function(meta, ...) UseMethod("get.data.matrix")
 
get.data.matrix.bayesTFR.mcmc.meta <- function(meta, ...) return (meta$tfr_matrix)

"get.countries.index" <- function(meta, ...) UseMethod("get.countries.index")

get.countries.index.bayesTFR.mcmc.meta  <- function(meta, ...) 
	return (if(is.null(meta$phase) || (meta$phase==2)) meta$id_DL else meta$id_phase3)

"get.countries.table" <- function(object, ...) UseMethod("get.countries.table")
get.countries.table.bayesTFR.mcmc.set <- function(object, ...) {
	ctable <- data.frame(code=object$meta$regions$country_code, name=object$meta$regions$country_name)
	if(!is.null(object$meta$phase) && (object$meta$phase==3)) ctable <- ctable[object$meta$id_phase3,]
	return(ctable)
}

get.countries.table.bayesTFR.prediction <- function(object, ...) {
	n <- dim(get.data.imputed(object))[2]
	return(data.frame(code=object$mcmc.set$meta$regions$country_code[1:n], 
					name=object$mcmc.set$meta$regions$country_name[1:n]))
}

