`startParallel` <- 
function(parallel.config, 
	process,
	qr.taus) {
	
	if (any(toupper(parallel.config[['BACKEND']]) == 'MULTICORE' | toupper(parallel.config[['BACKEND']]) == 'SNOW')) {
		stop(paste('\n\t', parallel.config[['BACKEND']], "no longer supported.  Please use the 'PARALLEL' package backend and R > 2.12 for parallel computation.\n"))
	}
	
	if (toupper(parallel.config[['BACKEND']]) == 'FOREACH' && (parallel.config[['TYPE']] != "doParallel" & !is.na(parallel.config[['TYPE']]))) {
			stop(paste('\n\t', parallel.config[['TYPE']], "no longer supported.  Please use doParallel and R > 2.12 for parallel computation.\n"))
	}
	
	workers <- NULL; par.type <- 'OTHER'; TAUS.LIST <- NULL

	if (!is.null(parallel.config[['CLUSTER.OBJECT']])) {
		if (!missing(qr.taus)) {
			workers <- length(eval(parse(text=parallel.config[['CLUSTER.OBJECT']])))
			chunk.size <- ceiling(length(qr.taus) / workers)
			TAUS.LIST <- vector("list", workers)
			for (chunk in 0:(workers-1)) {
				lower.index <- chunk*chunk.size+1
				upper.index <- min((chunk+1)*chunk.size, length(qr.taus))
				TAUS.LIST[[chunk+1]] <- qr.taus[lower.index:upper.index]
			}
		}
		clusterEvalQ(eval(parse(text=parallel.config[['CLUSTER.OBJECT']])), library(SGP))
		par.start <- list(internal.cl=eval(parse(text=parallel.config[['CLUSTER.OBJECT']])), par.type='SNOW')
		clusterExport(eval(parse(text=parallel.config[['CLUSTER.OBJECT']])), "par.start", envir=2)
		return(list(internal.cl=eval(parse(text=parallel.config[['CLUSTER.OBJECT']])), 
			par.type='SNOW', TAUS.LIST=TAUS.LIST))
	}
	
	###  Basic checks - default to ANY percentiles or projections WORKERS.
	
	if (is.numeric(parallel.config[['WORKERS']])) {
		message(paste("\t\tNOTE: ", process, " workers not specified.  Numeric value from WORKERS (", parallel.config[['WORKERS']], ") will be used for all processes.\n", sep=""))
		parallel.config[['WORKERS']][[process]] <- parallel.config[['WORKERS']]
	}
	if (is.null(parallel.config[['WORKERS']][[process]])) {
		if (!is.null(parallel.config[['WORKERS']])) {
			 tmp.indx <- grep(strsplit(process, "_")[[1]][2], names(parallel.config[['WORKERS']]))
			 if (any(!is.na(tmp.indx))) {
				 parallel.config[['WORKERS']][[process]] <- parallel.config[['WORKERS']][[tmp.indx]]
				 message(paste("\t\tNOTE: ", process, "workers not defined specifically.", names(parallel.config[['WORKERS']][tmp.indx]), 
				 	"WORKERS will be used  (", parallel.config[['WORKERS']][tmp.indx], "worker processors)."))
			 }
		} # See if still NULL and stop:
		if (is.null(parallel.config[['WORKERS']][[process]])) {
			# stop(paste(process, "workers must be specified."))
			parallel.config[['WORKERS']][[process]] <- 1
			message("\n\t\tNOTE: ", process, " workers not specified!  WORKERS will be set to a single (1) process.\n", sep="")
		}
	}
	
	Lower_Level_Parallel <- NULL
	if (all(c("PERCENTILES", "TAUS") %in% names(parallel.config[['WORKERS']]))) {
		# if (as.numeric(parallel.config[['WORKERS']][['PERCENTILES']])==1) {
		# 	Lower_Level_Parallel <- parallel.config
		# } else stop("Both TAUS and PERCENTILES can not be executed in Parallel at the same time.")
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both TAUS and PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		Lower_Level_Parallel <- parallel.config
	}
	if (all(c("PERCENTILES", "SIMEX") %in% names(parallel.config[['WORKERS']]))) {
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both SIMEX and PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		Lower_Level_Parallel <- parallel.config
	}
	
	if (all(c("BASELINE_PERCENTILES", "TAUS") %in% names(parallel.config[['WORKERS']]))) {
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both TAUS and BASELINE_PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		Lower_Level_Parallel <- parallel.config
	}
	if (all(c("BASELINE_PERCENTILES", "SIMEX") %in% names(parallel.config[['WORKERS']]))) {
		if (.Platform$OS.type != "unix" | "SNOW_TEST" %in% names(parallel.config)) stop("Both SIMEX and BASELINE_PERCENTILES can not be executed in Parallel at the same time in Windows OS or using SNOW type backends.")
		Lower_Level_Parallel <- parallel.config
	}

	###  Basic configuration
	
	if (toupper(parallel.config[['BACKEND']]) == 'FOREACH') {
		if (!is.na(parallel.config[['TYPE']]) & !identical(parallel.config[['TYPE']], "NA")) {
			eval(parse(text=paste("suppressPackageStartupMessages(require(", parallel.config[['TYPE']], "))")))
		} else parallel.config[['TYPE']] <- "doParallel"

		# if (parallel.config[['TYPE']]=="doMC" & is.null(parallel.config[['OPTIONS']][["preschedule"]])) {
			# if (is.list(parallel.config[['OPTIONS']])) {
				# parallel.config[['OPTIONS']][["preschedule"]]=FALSE
			# }	else parallel.config[['OPTIONS']]=list(preschedule=FALSE)
		# }

		if (parallel.config[['TYPE']]=="doParallel") { 
			if (.Platform$OS.type == "unix" & par.type == "OTHER") par.type <- 'MULTICORE' 
			if (.Platform$OS.type != "unix" & par.type == "OTHER") par.type <- 'SNOW'
			if (par.type == 'MULTICORE' & is.null(parallel.config[['OPTIONS']][["preschedule"]])) {
				if (is.list(parallel.config[['OPTIONS']])) {
					parallel.config[['OPTIONS']][["preschedule"]]=FALSE
				}	else parallel.config[['OPTIONS']]=list(preschedule=FALSE)
			}
		} # END doParallel
		
		foreach.options <- parallel.config[['OPTIONS']] # works fine if NULL
	} #  END FOREACH

	# if (toupper(parallel.config[['BACKEND']]) == 'MULTICORE') {
		# par.type <- 'MULTICORE'
	# }

	# if (toupper(parallel.config[['BACKEND']]) == 'SNOW') {
		# par.type <- 'SNOW'
	# }

	if (toupper(parallel.config[['BACKEND']]) == 'PARALLEL') {
		if (is.null(parallel.config[['TYPE']]) & !is.null(parallel.config[['SNOW_TEST']])) parallel.config[['TYPE']] <- 'PSOCK'
		if (!is.null(parallel.config[['TYPE']])) {
			if (!parallel.config[['TYPE']] %in% c('SOCK', 'PSOCK', 'MPI')) {
				stop("The 'snow' package will be used when 'parallel.config$TYPE' is specified and BACKEND=='PARALLEL'.  List element must be 'SOCK' ('PSOCK') or 'MPI'.")
			}
			par.type <- 'SNOW'
		} else {
			if (.Platform$OS.type == "unix") par.type <- 'MULTICORE' 
			if (.Platform$OS.type != "unix") par.type <- 'SNOW'; parallel.config[['TYPE']] <- 'PSOCK'
		}
	}
	
	if (par.type == 'SNOW') {
		if (is.null(parallel.config[['TYPE']])) stop("The 'parallel.config$TYPE' must be specified ('PSOCK' or 'MPI')")
		if (!parallel.config[['TYPE']] %in% c('PSOCK','MPI', 'doParallel')) stop("The 'parallel.config$TYPE' must be 'PSOCK', 'MPI' or 'doParallel'.")
	}


	###  Set up workers and spin up clusters / register workers
	
	if (!is.null(parallel.config[['WORKERS']][[process]])) {
		workers <- parallel.config[['WORKERS']][[process]]
	} else workers <- parallel.config[['WORKERS']]
	if (is.null(workers)) workers <- getOption("cores")
	if (is.null(workers)) stop("parallel.config$WORKERS must, at a minimum, contain the number of parallel workers for all processes, 
		or getOption('cores') must be specified to use MULTICORE parallel processing.")

	###
	###  Need this for all flavors - move to startParallel
	###

	if (process=='TAUS') {
		if (workers > 3) {
			if (workers %in% 4:10) {
				tmp.sml <- ceiling((length(qr.taus) / workers)*0.75)
				tmp.lrg <- ceiling((length(qr.taus)-(2*tmp.sml))/(workers-2))
				chunk.size <- c(tmp.sml, rep(tmp.lrg, (workers-2)), tmp.sml)
				if (sum(chunk.size) > length(qr.taus)) {
					over <- (sum(chunk.size) - length(qr.taus)); index <- 0
					while(over != 0) {
						if (over %% 2 == 0) {
							index <- index + 1
							chunk.size[(length(chunk.size)-(index))] <- chunk.size[(length(chunk.size)-(index))]-1
						} else chunk.size[(index + 1)] <- chunk.size[(index + 1)]-1
						over <- over - 1
					}
				}
			}
			if (workers > 10) {
				tmp.sml.a <- ceiling((length(qr.taus) / workers)*0.334)
				tmp.sml.b <- ceiling((length(qr.taus) / workers)*0.666)
				tmp.lrg <- ceiling((length(qr.taus)-(2*sum(tmp.sml.a, tmp.sml.b)))/(workers-4))
				chunk.size <- c(tmp.sml.a, tmp.sml.b, rep(tmp.lrg, (workers-4)), tmp.sml.b, tmp.sml.a)
				if (sum(chunk.size) > length(qr.taus)) {
					over <- (sum(chunk.size) - length(qr.taus)); index <- 0
					while(over != 0) {
						if (over %% 2 != 0) {
							index <- index + 1
							chunk.size[(length(chunk.size)-(index + 1))] <- chunk.size[(length(chunk.size)-(index + 1))]-1
						} else chunk.size[(index + 2)] <- chunk.size[(index + 2)]-1
						over <- over -1
					}
				}
			}
		}	else chunk.size <- rep(ceiling(length(qr.taus) / workers), workers)

		TAUS.LIST <- vector("list", workers)
		count <- index <- 1
		for (ch in chunk.size) {
			TAUS.LIST[[index]] <- qr.taus[count:(count+ch-1)]
			count <- (count+ch); index <- index + 1
		}
		if (sum(chunk.size) > length(qr.taus)) for(l in 1:length(TAUS.LIST))  TAUS.LIST[[l]] <- TAUS.LIST[[l]][!is.na(TAUS.LIST[[l]])]
	}
	
	###
	### END to startParallel
	###
	
	if (toupper(parallel.config[['BACKEND']]) == 'FOREACH') {
		if (parallel.config[['TYPE']]=="NA") {
			registerDoSEQ() # prevents warning message
			return(list(foreach.options=foreach.options, par.type=par.type, TAUS.LIST=TAUS.LIST))
		}
		# if (parallel.config[['TYPE']]=="doMC") {
			# registerDoMC(workers)
			# return(list(foreach.options=foreach.options, par.type=par.type, TAUS.LIST=TAUS.LIST))
		# }
		# if (parallel.config[['TYPE']]=='doMPI') {
			# doPar.cl <- startMPIcluster(count=workers)
			# registerDoMPI(doPar.cl)
			# return(list(doPar.cl=doPar.cl, foreach.options=foreach.options, par.type=par.type))
		# }
		# if (parallel.config[['TYPE']]=='doRedis') {
			# redisWorker('jobs', port=10187) #  Doesn't seem to work.  Maybe get rid of this option/flavor?
			# registerDoRedis('jobs')
			# startLocalWorkers(n=workers, queue='jobs')
			# return(list(jobs='jobs', foreach.options=foreach.options, par.type=par.type))
		# }
		# if (parallel.config[['TYPE']]=='doSNOW') {
			# doPar.cl=makeCluster(workers, type='SOCK')
			# registerDoSNOW(doPar.cl)
			# return(list(doPar.cl=doPar.cl, foreach.options=foreach.options, par.type=par.type))
		# }
		if (!is.null(parallel.config[['SNOW_TEST']])) par.type <- 'SNOW' # To test SNOW on Linux
		if (parallel.config[['TYPE']]=="doParallel") {
			if (par.type == 'SNOW') {
				doPar.cl <- makePSOCKcluster(workers)
				registerDoParallel(doPar.cl)
				clusterEvalQ(doPar.cl, library(SGP))
				# foreach.options <- list(attachExportEnv=TRUE)
				return(list(doPar.cl=doPar.cl, foreach.options=foreach.options, par.type=par.type, TAUS.LIST=TAUS.LIST))
			} else {
				registerDoParallel(workers)
				return(list(foreach.options=foreach.options, par.type=par.type, TAUS.LIST=TAUS.LIST))
			}
		}
	} # END if (FOREACH)

	if (par.type=='SNOW') {
		# if (parallel.config[['TYPE']]=='MPI') {
			# if (exists('par.start')) return() #don't try to restart a new config
		# }
		internal.cl <- makeCluster(eval(parse(text=workers)), type=parallel.config[['TYPE']]) # eval workers in case 'names' used
		clusterEvalQ(internal.cl, library(SGP))
		return(list(internal.cl=internal.cl, par.type=par.type, TAUS.LIST=TAUS.LIST)) #  workers=workers,
	}

	if (par.type=='MULTICORE') {
		return(list(workers=workers, par.type=par.type, TAUS.LIST=TAUS.LIST, Lower_Level_Parallel=Lower_Level_Parallel))
	}
} ### END startParallel Funtion
