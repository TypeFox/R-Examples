pcit <- function(m, force.serial=FALSE, force.parallel=FALSE, nslaves=NULL, verbose=getOption("verbose"), tol.type=c("mean", "min", "max", "median"), pass.type=c("file", "memory", "db") ) {
	pass.type <- match.arg(pass.type)
	tol.type <- match.arg(tol.type)
	
	# Create data structure to store the results
	results <- list()
	results$idx <- NULL
	
	if(force.serial && force.parallel) {
		cat("  WARN: You specified both force.serial=TRUE and force.parallel=TRUE. I'll assume you want me to choose which implementation to use!\n")
	} else if(force.parallel) {
		# attempt to setup parallel environment
		run_parallel <- .setup_mpi()
		if(!run_parallel) { stop("  Unable to force the parallel implementation of pcit()") }
	} else if(force.serial) {
		run_parallel <- FALSE
	} else {
		# neither force.serial or force.parallel set to true
		# attempt to autodetect parallel environment
		run_parallel <- .setup_mpi()
	}
	
	
	if (run_parallel) {
		if(verbose) { cat("  Running parallel implementation of pcit().\n") }
	} else {
		if(verbose) { cat("  Running serial implementation of pcit().\n") }
		results <- .pcit(m, verbose=verbose, tol.type=tol.type)
		results$idx <- .sub2ind(results$idx, nrow=nrow(m), ncol=ncol(m))
		return(results)
	}
	
	# we only get here if we have a detected a parallel environment
	# We're gonna pass the data to the slaves using one of the possible options
	switch(pass.type,
			file = {
				RData <- tempfile()
				save(m, file=RData)
			},
			memory = {
				cat("WARN: setting pass.file=FALSE may provide a small amount of speedup by passing data around in memory.\n  However, there is a limit to the size of the data set that can be passed in memory.\n  If you receive an 'Error: serialization is too large to store in a raw vector', you have reached this limit and should use pass.type=\"file\" (the default)\n")
			},
			db = {
				stop("Using a DB to pass data is not yet implemented.\n")
			}
	)
	
	nSlaves <- Rmpi::mpi.comm.size()-1
	
	# Function the slaves will call to perform a validation on the
	# fold equal to their slave number.
	# Assumes the following object(s) will be passed to it by the master: m
	slavefunction <- function() {
		# load the pcit library in each slave
		require('PCIT', quietly=TRUE)
		
		# Note the use of the tag for sent messages: 
		#     1=ready_for_task, 2=done_task, 3=exiting 
		# Note the use of the tag for received messages: 
		#     1=task, 2=done_tasks 
		junk <- 0
		done <- 0
		while (done != 1) {
			# Signal being ready to receive a new task 
			Rmpi::mpi.send.Robj(junk,0,1)
			
			# Receive a object list from master
			objList <- Rmpi::mpi.recv.Robj(Rmpi::mpi.any.source(),Rmpi::mpi.any.tag())
			task_info <- Rmpi::mpi.get.sourcetag()
			tag <- task_info[2]
			
			if (tag == 1) {
				# do a task, taking data/objects from the objList
				# Find thos entries that are zero'd out in the m_sub_pcit matrix and convert them back to indices for the full matrix
				# Lets not check to see see which zero's were changed only by this slave as we'd need more time and memory to do so - i.e. a nGenes x nGenes Boolean matrix for each check
				switch(pass.type,
						file = {
							# load the whole data matrix from file and run PCIT() on it
							load(RData)
							results <- .pcit(m, x=objList$x, verbose=verbose, tol.type=tol.type)
							results$nonsignifOffset <- 1
						},
						memory = {
							results <- .pcit(m, x=objList$x, verbose=verbose, tol.type=tol.type)
							results$nonsignifOffset <- 1
						},
						db = {
							stop("Using a DB to pass data is not yet implemented.\n")
						}
				)
				
				results$idx <- .sub2ind(results$idx, nrow=nrow(m), ncol=ncol(m))
				
				# Send a results message back to the master
				Rmpi::mpi.send.Robj(results,0,2)
				
			} else if (tag == 2) {
				# maybe send date() info back to the master so it can compile debugging info on the distribution of times each slave took?
				done <- 1
				if (verbose) {
					cat("Slave", Rmpi::mpi.comm.rank(), "- FINISHED:", date(), "\n")
				}
			}
			# We'll just ignore any unknown messages
		}
		
		Rmpi::mpi.send.Robj(junk,0,3)
	}
	
	# We're in the parent.  
	# send data object(s) to all the slaves
	Rmpi::mpi.bcast.Robj2slave(tol.type)
	
	switch(pass.type,
			file = {
				Rmpi::mpi.bcast.Robj2slave(RData)
			},
			memory = {
				# This means each slave has it's own copy of the master correlation matrix
				Rmpi::mpi.bcast.Robj2slave(m)
			},
			db = {
				stop("Using a DB to pass data is not yet implemented.\n")
			}
	)
	
	# Send function(s) to all the slaves
	Rmpi::mpi.bcast.Robj2slave(slavefunction)
	# Call the function in all the slaves to get them ready to
	# undertake tasks
	Rmpi::mpi.bcast.cmd(slavefunction())
	
	# define a list of tasks
	tasks <- defineTasks(n=nrow(m), nSlaves=nSlaves)
	
	junk <- 0
	exited <- 0
	child_errors <- FALSE
	while (exited < nSlaves) { 
		# Receive a message from a slave 
		message <- Rmpi::mpi.recv.Robj(Rmpi::mpi.any.source(),Rmpi::mpi.any.tag()) 
		message_info <- Rmpi::mpi.get.sourcetag() 
		slave_id <- message_info[1] 
		tag <- message_info[2] 
		
		if (tag == 1) {
			# slave is ready for a task. Give it the next task, or tell it tasks 
			# are done if there are none.
			# If errors have been found in children, don't give any new tasks
			if (child_errors) {
				cat("WARN: Telling slave ",slave_id," to exit\n");
				Rmpi::mpi.send.Robj(junk,slave_id,2)
			}
			if (length(tasks) > 0) {
				switch(pass.type,
						file = {
							# send the job with data loaded from file
							Rmpi::mpi.send.Robj(tasks[[1]], slave_id, 1)
						},
						memory = {
							Rmpi::mpi.send.Robj(tasks[[1]], slave_id, 1)
						},
						db = {
							stop("Using a DB to pass data is not yet implemented.\n")
						}
				)
				
				tasks[[1]] <- NULL 
				
			} else {
				Rmpi::mpi.send.Robj(junk, slave_id, 2)
			}
			
		} else if (tag == 2) { 
			# The message contains results
			# compile each subset into a master result
			
			# for linear ind which are positions of non-significnat connections - basically, a zero in one/both of results and message
			if(verbose) {
				cat("Master: calculating unique indicies ... ")
			}
			
			if (is.null(results$idx)) {
				# if the master doesn't have any indices yet, just assign the result from the slave
				# The slave will only have zero'd out connections within the sub matrix it was dealing with
				results$idx <- message$idx
			} else {
				# combining data from slave to that which the master already knows about
				results$idx <- intersect(results$idx, message$idx)
			}
			
			if(verbose) {
				cat("DONE\n")
			}
			
			# for linear ind which are positions of significant connections - basically, a 1 in both results and message
			# BUT - only from the same overlapping sub matrix - need to append those indices from outside the overlapping region of the larger matrix
			# the following DOES NOT work!
			#results$nonsignif <- c(results$nonsignif, message$nonsignif)[duplicated(c(results$nonsignif, message$nonsignif))]
			
			#results$nonsignif <- append(results$nonsignif, message$nonsignif)
			
		} else if (tag == 3) {
			# A slave has closed down.
			exited <- exited + 1
		} else if (tag == 4) {
			exited <- exited + 1
			child_errors <- TRUE
			cat("ERROR: Slave ", slave_id, ": ", message,"\n")
		}
	}
	if(verbose) {
		cat("Master: completed pcit().\n")
	}
	# clean up
	switch(pass.type,
			file = {
				unlink(RData)
			},
			memory = {
				
			},
			db = {
				stop("Using a DB to pass data is not yet implemented.\n")
			}
	)
	if (child_errors) {
		stop("Errors in R slaves.  See last.warnings or check slave logs")
	}
	
	.freeSlaves()
	
	return(results)
}
