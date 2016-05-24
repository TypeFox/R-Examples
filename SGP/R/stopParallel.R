`stopParallel` <- 
function(parallel.config,
	par.start) {
	
	if (toupper(parallel.config[['BACKEND']]) == 'FOREACH') {
		if (is.na(parallel.config[['TYPE']])) parallel.config[['TYPE']] <- "NA" # Stop checks on NA below.

		# if (parallel.config[['TYPE']]=="doMPI") {
			# closeCluster(par.start$doPar.cl)
			# return()
		# }

		if (parallel.config[['TYPE']]=="doParallel" & par.start$par.type=='SNOW') {
			stopCluster(par.start$doPar.cl)
			return()
		}
		
		# if (parallel.config[['TYPE']]=="doSNOW") {
			# stopCluster(par.start$doPar.cl)
			# return()
		# }
		
		# if (parallel.config[['TYPE']]=="doRedis") {
			# removeQueue(par.start$jobs)
			# return()
		# }
	} #  END FOREACH

	# Nothing required for MULTICORE (or doMC)
	# if (toupper(parallel.config[['BACKEND']]) == 'MULTICORE') {
	# }

	if (par.start$par.type == 'SNOW') {
		if (is.null(parallel.config[['CLUSTER.OBJECT']]))	 {
			stopCluster(par.start$internal.cl)
		}
	}
} ### END stopParallel Function
