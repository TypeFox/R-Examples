.setup_mpi <- function(nslaves=NULL) {
	# Default to running serially.  Only run in parallel if Rmpi
	# is installed AND there's multiple slaves AND pcit is called
	# from the master.
	run_parallel <- FALSE
	
	# Load the MPI Environment if not already there.
	if (!is.loaded('mpi_initialize')) {
		if (sum((.packages())=='Rmpi') == 0) {
			return(run_parallel)
		}
		requireNamespace('Rmpi', quietly=TRUE)
	}
	
	# Now, if Rmpi is loaded, make sure we have a bunch of slaves
	if (is.loaded('mpi_initialize')) {
		if (Rmpi::mpi.comm.size() > 2) {
			# It's already loaded and set up
			return(TRUE)
		}
		
		# Spawn as many slaves as possible
		if (Rmpi::mpi.comm.size() < 2) {
			if(is.null(nslaves)) {
				answer <- try(Rmpi::mpi.spawn.Rslaves(), silent=TRUE)
			} else {
				answer <- try(Rmpi::mpi.spawn.Rslaves(nslaves=nslaves), silent=TRUE)
			}
			if (class(answer) == "try-error") {
				cat("WARNING: Managed to load the Rmpi library but failed to spawn slaves (error message below) - Falling back to serial implementation.\n\n", geterrmessage())
				return(FALSE)
			}
		}
		
		# If there's multiple slaves, AND this process is the master,
		# run in parallel.
		if (Rmpi::mpi.comm.size() > 2) {
			if (Rmpi::mpi.comm.rank() == 0) {
				run_parallel <- TRUE
			}
		} else {
			cat("WARN: I could only find", Rmpi::mpi.comm.size()-1, "slaves to work with. Parallel programming usually requires >= 2 slaves, but I'll continue anyway!\n")
			run_parallel <- TRUE
		}
	}
	
	return(run_parallel)
}

.pcit <- function(m, x=1:(nrow(m)-2), tol.type=c("mean", "min", "max", "median"), verbose=getOption("verbose")) {
	# m is a matrix of direct correlations between genes using:
	# m <- similarity(d, type="none")
	
	if (! isSymmetric(m) ) {
		stop("The pcit() function requires a square, symmetrical matrix as input")
	}
	
	# use "mean" as the default for calculating the tolerances
	tol.type <- match.arg(tol.type)
	switch(tol.type,
		"mean" = { tol.type <- 1 },
		"min" = { tol.type <- 2 },
		"max" = { tol.type <- 3 },
		"median" = {
			tol.type <- 4
			stop("Tolerance calculation using median is not implemented in the PCIT Fortran code.")
		}
	)
	# set NA's to zero
	index <- is.na(m)
	if (sum(index)>0) {
		m[index] <- 0
		warning(paste(sum(index), " values were found to be N/A and set to zero before running pcit().", sep=""))
	}
	
	m_partials <- m
	
	result <- .Fortran( "pcit", correlations=as.single(m), partial_correlations=as.single(m_partials), nGenes=as.integer(nrow(m)), xVals=as.integer(x), nXVals=as.integer(length(x)), tolType=as.integer(tol.type) )
	
	result$partial_correlations <- matrix(result$partial_correlations,nrow(m),ncol(m))
	result$idx <- which(result$partial_correlations != 0, arr.ind=TRUE)
	
	# remove items from the results list which we don't want to return to the calling function
	result[c("correlations", "partial_correlations", "nGenes", "xVals", "nXVals", "tolType")] <- NULL
	
	# return only the idecies for those that PCIT found to be significant
	return(result)
}

.sub2ind <- function(x, y, nrow, ncol=NULL) {
	## Returns a linear index for the (x,y) coordinates passed in.
	if (is.matrix(x) || is.data.frame(x)) {
		stopifnot(ncol(x) == 2)
		if (!missing(y)) {
			if (missing(nrow)) {
				nrow <- y
			} else {
				ncol <- nrow
				nrow <- y
			}
		}
		y <- x[,2]
		x <- x[,1]
	}
	
	if (is.matrix(nrow)) {
		d <- dim(nrow)
		nrow <- d[1]
		ncol <- d[2]
	} else if (is.null(ncol)) {
		stop("Dimensions of matrix under-specified")
	}
	
	# Sanity check to ensure we got each var doing what it should be doing
	if (length(x) != length(y) || length(nrow) != 1 || length(ncol) != 1) {
		stop("I'm confused")
	}
	
	((x - 1) + ((y - 1) * nrow)) + 1
	
}

.freeSlaves <- function(mpi.exit=FALSE) {
	if (is.loaded("mpi_initialize")){
		if (Rmpi::mpi.comm.size(1) > 0){
			Rmpi::mpi.close.Rslaves()
		}
		if(Rmpi::mpi.exit) {
			Rmpi::mpi.exit()
		}
	}
}

.Last <- function(){
	.freeSlaves(mpi.exit=TRUE)
}
