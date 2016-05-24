if(getRversion() >= "2.15.1") utils::globalVariables(c("counter", "counter3"))

store.mcmc <- local({
	# Writes parameter values into ascii files - one file per parameter and country (if country-specific)
	##########################
	par.names <- c(tfr.parameter.names(trans=FALSE))
	par.cs.names <- c(tfr.parameter.names.cs(trans=FALSE, back.trans=FALSE), 'eps_T')
	var.names <- list(gamma='gamma_ci', d='d_c', Triangle_c4='Triangle_c4', eps_T='eps_Tc', U='U_c')
	
	default.buffer.size <- 10
	buffer <- buffer.cs <- NULL
	
	get.gamma <- function(mcmc, country) {
		return(mcmc$gamma_ci[country,])
	}
    get.eps_T <- function(mcmc, country) {
        return(t(mcmc$eps_Tc[,country]))
    }
	special.case <- c('gamma', 'eps_T')
	
	buffers.insert <- function(mcmc, countries=NULL) {
		counter <<- counter + 1
		if (is.null(countries)) {
			for (par in par.names) {
				if (is.element(par, mcmc$dontsave)) next
				buffer[[par]][counter,] <<- mcmc[[par]]
			}
			country.index <- mcmc$meta$id_DL
		} else {
			country.index <- countries
		}
		for (par in par.cs.names) {
			if (is.element(var.names[[par]], mcmc$dontsave)) next
			
			for (country in country.index){
				if (is.element(par, special.case)) {
					result <- eval(call(paste('get', par, sep='.'), mcmc, country))
				} else {
					result <- mcmc[[var.names[[par]]]][country]
				}
				buffer.cs[[par]][[country]][counter,] <<- result
			}
		}
	}
		
	buffers.ini <- function(mcmc, size, countries=NULL) {
		buffer <<- list()
		if (is.null(countries)) {
			for (par in par.names) {
				if (is.element(par, mcmc$dontsave)) next
				buffer[[par]] <<- matrix(NA, ncol=length(mcmc[[par]]), nrow=size)
			}
			country.index <- mcmc$meta$id_DL
		} else {
			country.index <- countries
		}
		buffer.cs <<-list()
		for (par in par.cs.names) {
			if (is.element(var.names[[par]], mcmc$dontsave)) next
			buffer.cs[[par]] <<- list()
			for (country in country.index){
				if (is.element(par, special.case)) {
					v <- eval(call(paste('get', par, sep='.'), mcmc, country))
				} else {
					v <- mcmc[[var.names[[par]]]][country]
				}
				buffer.cs[[par]][[country]] <<- matrix(NA, ncol=length(v), nrow=size)
			}
		}
		counter <<- 0
	}
	
	
	do.flush.buffers <- function(mcmc, append=FALSE, countries=NULL, verbose=FALSE) {
		if (verbose)
			cat("Flushing results into disk.\n")
		output.dir <- file.path(mcmc$meta$output.dir, mcmc$output.dir)
		if(!file.exists(output.dir)) 
			dir.create(output.dir)
		open <- if(append) 'a' else 'w'
		if (is.null(countries)) {
			for(par in par.names) { # write country-independent parameters
				if (is.null(buffer[[par]])) next
				if (counter == 1) {
					values <- t(buffer[[par]][1:counter,])
				} else {
					values <- buffer[[par]][1:counter,]
				}
				write.values.into.file.cindep(par, values, output.dir, mode=open, 
												compression.type=mcmc$compression.type)
			}
			country.index <- mcmc$meta$id_DL	
		} else {
			country.index <- countries
		}
		for (par in par.cs.names) { # write country-specific parameters
			if (is.null(buffer.cs[[par]])) next
			for (country in country.index){
				if (counter == 1) {
					values <- t(buffer.cs[[par]][[country]][1:counter,])
				} else {
					values <- buffer.cs[[par]][[country]][1:counter,]
				}
				write.values.into.file.cdep(par, values, output.dir, 
						get.country.object(country, meta=mcmc$meta, index=TRUE)$code, mode=open, 
											compression.type=mcmc$compression.type)
			}
		}
		resmc <- as.list(mcmc)
		class(resmc) <- 'bayesTFR.mcmc'
		store.bayesTFR.object(resmc, output.dir)
	}
	
	store <- function(mcmc, append=FALSE, flush.buffer=FALSE, countries=NULL, verbose=FALSE) {
		# If countries is not NULL, only country-specific parameters 
		# for those countries (given as index) are stored
		buffer.size <- mcmc$meta$buffer.size
		if (is.null(buffer.size)) buffer.size <- default.buffer.size
		if (is.null(buffer)) buffers.ini(mcmc, buffer.size, countries=countries)
		buffers.insert(mcmc, countries=countries)
		flushed <- FALSE
		if (flush.buffer || (counter >= buffer.size)) {
			do.flush.buffers(mcmc, append=append, countries=countries, verbose=verbose)
			buffer <<- buffer.cs <<- NULL
			flushed <- TRUE
		}
		return(flushed)
	}

})

store.mcmc3 <- local({
	# Writes parameter values into ascii files - one file per parameter and country (if country-specific)
	##########################
	par.names <- tfr3.parameter.names()
	par.cs.names <- tfr3.parameter.names.cs()
	
	default.buffer.size <- 10
	buffer3 <- buffer3.cs <- NULL
		
	buffers.insert <- function(mcmc, countries=NULL) {
		counter3 <<- counter3 + 1
		if (is.null(countries)) {
			for (par in par.names) buffer3[[par]][counter3,] <<- mcmc[[par]]
			country.index <- 1: mcmc$meta$nr.countries
		} else country.index <- countries
		for (par in par.cs.names) {			
			for (country in country.index)
				buffer3.cs[[par]][[country]][counter3,] <<- if(is.null(dim(mcmc[[par]]))) mcmc[[par]][country] 
                								          else mcmc[[par]][,country]
		}
	}
		
	buffers.ini <- function(mcmc, size, countries=NULL) {
		buffer3 <<- list()
		if (is.null(countries)) {
			for (par in par.names) 
				buffer3[[par]] <<- matrix(NA, ncol=length(mcmc[[par]]), nrow=size)
			country.index <- 1:mcmc$meta$nr.countries
		} else country.index <- countries
		buffer3.cs <<-list()
		for (par in par.cs.names) {
			buffer3.cs[[par]] <<- list()
			for (country in country.index){
				v <- if(is.null(dim(mcmc[[par]]))) mcmc[[par]][country] else mcmc[[par]][,country]
				buffer3.cs[[par]][[country]] <<- matrix(NA, ncol=length(v), nrow=size)
			}
		}
		counter3 <<- 0
	}
	
	do.flush.buffers <- function(mcmc, append=FALSE, countries=NULL, verbose=FALSE) {
		if (verbose)
			cat("Flushing results into disk.\n")
		output.dir <- file.path(mcmc$meta$output.dir, mcmc$output.dir)
		if(!file.exists(output.dir)) 
			dir.create(output.dir)
		open <- if(append) 'a' else 'w'
		if (is.null(countries)) {
			for(par in par.names) { # write country-independent parameters
				if (is.null(buffer3[[par]])) next
				values <- if (counter3 == 1) t(buffer3[[par]][1:counter3,])
				 			else buffer3[[par]][1:counter3,]
				write.values.into.file.cindep(par, values, output.dir, mode=open, 
												compression.type=mcmc$compression.type)
			}
			country.index <- 1:mcmc$meta$nr.countries	
		} else country.index <- countries

		for (par in par.cs.names) { # write country-specific parameters
			if (is.null(buffer3.cs[[par]])) next
			for (country in country.index){
				values <- if (counter3 == 1) t(buffer3.cs[[par]][[country]][1:counter3,])
							else buffer3.cs[[par]][[country]][1:counter3,]
				write.values.into.file.cdep(par, values, output.dir, 
						get.country.object(mcmc$meta$id_phase3[country], meta=mcmc$meta$parent, index=TRUE)$code, mode=open, 
											compression.type=mcmc$compression.type)
			}
		}
		resmc <- as.list(mcmc)
		class(resmc) <- 'bayesTFR.mcmc'
		store.bayesTFR.object(resmc, output.dir)
	}
	
	store <- function(mcmc, append=FALSE, flush.buffer=FALSE, countries=NULL, verbose=FALSE) {
		# If countries is not NULL, only country-specific parameters 
		# for those countries (given as index) are stored
		buffer.size <- mcmc$meta$buffer.size
		if (is.null(buffer.size)) buffer.size <- default.buffer.size
		if (is.null(buffer3)) buffers.ini(mcmc, buffer.size, countries=countries)
		buffers.insert(mcmc, countries=countries)
		flushed <- FALSE
		if (flush.buffer || (counter3 >= buffer.size)) {
			do.flush.buffers(mcmc, append=append, countries=countries, verbose=verbose)
			buffer3 <<- buffer3.cs <<- NULL
			flushed <- TRUE
		}
		return(flushed)
	}

})


.get.compression.settings <- function(compression.type='None') {
	if(is.null(compression.type)) compression.type <- 'None'
	return(switch(compression.type,
							None=c('file', '', ''),
							xz = c('xzfile', '.xz', 'b'),
							bz = c('bzfile', '.bz2','b'),
							gz = c('gzfile', '.gz', 'b')))
}

do.write.values.into.file <- function(filename, data, mode, compression.type='None') {
	cmd.suffix.mode <- .get.compression.settings(compression.type)
	#con <- bzfile(filename, open=mode)
	con <- do.call(cmd.suffix.mode[1], list(paste(filename, cmd.suffix.mode[2], sep=''), 
				open=paste(mode, cmd.suffix.mode[3], sep='')))
	write.table(data, file=con, row.names=FALSE, col.names = FALSE, sep=" ")
	close(con)
}

write.values.into.file.cindep <- function(par, data, output.dir, mode='w', compression.type='None') {
	do.write.values.into.file(file.path(output.dir, paste(par,'txt', sep='.')), data, mode=mode, 
									compression.type=compression.type)
}

write.table.into.file.cindep <- function(data, ...) {
	for (par in colnames(data))
		write.values.into.file.cindep(par, data[,par], mode='w', ...)
}

write.values.into.file.cdep <- function(par, data, output.dir, country.code, mode='w', compression.type='None') {
	do.write.values.into.file(file.path(output.dir, paste(par,"_country", country.code, ".txt",sep = "")), 
									data, mode=mode, compression.type=compression.type)
}

write.table.into.file.cdep <- function(data, ...) {
	for (par in colnames(data))
		write.values.into.file.cdep(par, data[,par], mode='w', ...)
}

store.bayesTFR.object <- function(mcmc, output.dir) {
	bayesTFR.mcmc <- mcmc
	for (item in bayesTFR.mcmc$dontsave)  # don't save meta and some other data
		bayesTFR.mcmc[[item]] <- NULL
	bayesTFR.mcmc$meta <- NULL
	save(bayesTFR.mcmc, file=file.path(output.dir, 'bayesTFR.mcmc.rda'))
}

store.bayesTFR.meta.object <- function(meta, output.dir) {
	bayesTFR.mcmc.meta <- meta
	save(bayesTFR.mcmc.meta, file=file.path(output.dir, 'bayesTFR.mcmc.meta.rda'))
}

store.bayesTFR.prediction <- function(pred, output.dir=NULL) {
	bayesTFR.prediction <- pred
	if (is.null(output.dir)) output.dir <- pred$output.directory
	save(bayesTFR.prediction, file=file.path(output.dir, 'prediction.rda'))
}

store.bayesTFR.convergence <- function(diag, thin, burnin, output.dir){
	save.file <- file.path(output.dir, paste('bayesTFR.convergence_', thin, '_', burnin, '.rda', sep=''))
	bayesTFR.convergence <- diag
	save(bayesTFR.convergence, file=save.file)
	return(save.file)
}