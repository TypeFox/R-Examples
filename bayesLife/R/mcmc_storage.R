if(getRversion() >= "2.15.1") utils::globalVariables("counter")

store.e0.mcmc <- local({
	# Writes parameter values into ascii files - one file per parameter and country (if country-specific)
    ##########################
    par.names <- e0.parameter.names()
    par.cs.names <- e0.parameter.names.cs()
        
    default.buffer.size <- 10
    buffer <- buffer.cs <- NULL
                
    buffers.insert <- function(mcmc, countries=NULL) {
		counter <<- counter + 1
		if (is.null(countries)) {
        	for (par in par.names) buffer[[par]][counter,] <<- mcmc[[par]]
			country.index <- 1: mcmc$meta$nr.countries
        } else country.index <- countries
        for (par in par.cs.names) {
			for (country in country.index)
                buffer.cs[[par]][[country]][counter,] <<- if(is.null(dim(mcmc[[par]]))) mcmc[[par]][country] 
                								          else mcmc[[par]][,country]
        }
	}
                
	buffers.ini <- function(mcmc, size, countries=NULL) {
		buffer <<- list()
		if (is.null(countries)) {
        	for (par in par.names) buffer[[par]] <<- matrix(NA, ncol=length(mcmc[[par]]), nrow=size)
			country.index <- 1: mcmc$meta$nr.countries
        } else country.index <- countries
		buffer.cs <<-list()
        for (par in par.cs.names) {
			buffer.cs[[par]] <<- list()
			for (country in country.index){
				v <- if(is.null(dim(mcmc[[par]]))) mcmc[[par]][country] else mcmc[[par]][,country]
                buffer.cs[[par]][[country]] <<- matrix(NA, ncol=length(v), nrow=size)
            }
        }
        counter <<- 0
	}
                
	do.flush.buffers <- function(mcmc, append=FALSE, countries=NULL, verbose=FALSE) {
		if (verbose)
			cat("Flushing results into disk.\n")
		output.dir <- file.path(mcmc$meta$output.dir, mcmc$output.dir)
		if(!file.exists(output.dir)) dir.create(output.dir)
		open <- if(append) 'a' else 'w'
		if (is.null(countries)) {
			for(par in par.names) { # write country-independent parameters
				if (is.null(buffer[[par]])) next
                values <- if (counter == 1) t(buffer[[par]][1:counter,]) else buffer[[par]][1:counter,]
                bayesTFR:::write.values.into.file.cindep(par, values, output.dir, mode=open, 
												compression.type=mcmc$compression.type)
            }
            country.index <- 1: mcmc$meta$nr.countries        
		} else country.index <- countries
		for (par in par.cs.names) { # write country-specific parameters
			if (is.null(buffer.cs[[par]])) next
            for (country in country.index){
            	values <- if (counter == 1) t(buffer.cs[[par]][[country]][1:counter,]) 
            				else values <- buffer.cs[[par]][[country]][1:counter,]
            parname <- par
			bayesTFR:::write.values.into.file.cdep(parname, values, output.dir, 
            		get.country.object(country, meta=mcmc$meta, index=TRUE)$code, mode=open, 
											compression.type=mcmc$compression.type)
            }
        }
        resmc <- as.list(mcmc)
		class(resmc) <- 'bayesLife.mcmc'
		store.bayesLife.object(resmc, output.dir)
	}
        
	store <- function(mcmc, append=FALSE, flush.buffer=FALSE, countries=NULL, verbose=FALSE) {
		# If countries is not NULL, only country-specific parameters 
		# for those countries (given as index) are stored
		buffer.size <- mcmc$meta$buffer.size
		if (is.null(buffer.size)) buffer.size <- default.buffer.size
		if (is.null(buffer)) buffers.ini(mcmc, buffer.size, countries=countries)
		buffers.insert(mcmc, countries=countries)
		if (flush.buffer || (counter >= buffer.size)) {
			do.flush.buffers(mcmc, append=append, countries=countries, verbose=verbose)
			buffer <<- buffer.cs <<- NULL
        }
	}

})

store.bayesLife.meta.object <- function(meta, output.dir) {
        bayesLife.mcmc.meta <- meta
        save(bayesLife.mcmc.meta, file=file.path(output.dir, 'bayesLife.mcmc.meta.rda'))
}

store.bayesLife.object <- function(mcmc, output.dir) {
        bayesLife.mcmc <- mcmc
        bayesLife.mcmc$meta <- NULL
        save(bayesLife.mcmc, file=file.path(output.dir, 'bayesLife.mcmc.rda'))
}

store.bayesLife.prediction <- function(pred, output.dir=NULL) {
	bayesLife.prediction <- pred
	if (is.null(output.dir)) output.dir <- pred$output.directory
	save(bayesLife.prediction, file=file.path(output.dir, 'prediction.rda'))
}

store.bayesLife.convergence <- function(diag, thin, burnin, output.dir){
	save.file <- file.path(output.dir, paste('bayesLife.convergence_', thin, '_', burnin, '.rda', sep=''))
	bayesLife.convergence <- diag
	save(bayesLife.convergence, file=save.file)
	return(save.file)
}