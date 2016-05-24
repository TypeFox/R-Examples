runjags.readin <- function(directory, silent.jags, target.adapt, target.burnin, target.iters, n.chains, monitor, method, method.options, suspended=FALSE, showoutput=FALSE, read.monitor=NA, sub.samples=FALSE, sub.chains=FALSE, force.read=FALSE){
	
	### force.read is currently ignored apart from to give a warning but may be implemented later

	save.directory <- getwd()
	temp.directory <- directory
	
	unfinished=allok=forcefailed <- FALSE
	on.exit({

		if (!allok || forcefailed){
			
			failedjagsmodel <- paste(readLines("model.txt", warn=FALSE), collapse="\n")
			class(failedjagsmodel) <- "runjagsmodel"
			assign("model", failedjagsmodel, envir=failedjags)

			failedd <- paste(readLines("data.txt", warn=FALSE), collapse="\n")
			class(failedd) <- "runjagsdata"
			assign("data", failedd, envir=failedjags)
			
			failedi <- character(orig.n.chains)
			for(i in 1:orig.n.chains){
				failedi[i] <- paste(readLines(paste("inits",i,".txt",sep=""), warn=FALSE), collapse="\n")
			}
			class(failedi) <- "runjagsinits"
			assign("inits", failedi, envir=failedjags)
			
		}
		
		if(is.na(showoutput)){
			if(suspended && unfinished) showoutput <- TRUE else showoutput <- FALSE
		}
		# We also want to save and print jagsoutput if it's a background JAGS object, or if we want to showoutput:
		if(!allok || suspended || showoutput || forcefailed){
			failedo <- character(orig.n.sims)
			for(i in 1:orig.n.sims){
				failedo[i] <- paste(readLines(paste("sim.",i,"/jagsoutput.txt",sep=""), warn=FALSE), collapse="\n")
			}
			names(failedo) <- paste('Process',1:orig.n.sims,sep='_')
			class(failedo) <- "runjagsoutput"
			if(showoutput) print(failedo)
			assign("output",failedo, envir=failedjags)
		}
		
		setwd(save.directory)
	})
	
	setwd(temp.directory)

	orig.n.chains <- n.chains

	# dummy variable to get rid of binding warnings - actually loaded from simchainsinfo:
	sim.chains <- list()
	success <- try(load("simchainsinfo.Rsave"))
	if(inherits(success, 'try-error')) stop("The required 'simchainsinfo.Rsave' file was not found in the root simulation directory, please file a bug report to the package developer!")
	orig.n.sims=n.sims <- length(sim.chains)
	orig.sim.chains <- sim.chains
	
  	if(length(list.files(pattern="JAGS.out",recursive=TRUE))>0){
   		return(list(errormessage="You are using a version of JAGS prior to 0.99.0, which is no longer supported.  Please update JAGS and try again"))
	}
	
	# Sort out behaviour for crashed sims:
	part.return <- FALSE
	which.chains <- sub.chains
	if(identical(sub.chains, TRUE)){
		which.chains <- 1:n.chains
		part.return <- TRUE
	}
	if(identical(sub.chains, FALSE) || length(sub.chains)==0){
		which.chains <- 1:n.chains
		part.return <- FALSE
	}
	if(any(which.chains > n.chains)){
		warning('The which.chains argument specified more than the available number of chains', call.=FALSE)
		which.chains <- which.chains[which.chains <= n.chains]
	}
	if(n.chains!=length(which.chains))
		part.return <- TRUE
	n.chains <- length(which.chains)
	
	simfolders <- paste("sim.",1:n.sims,sep="")

	if(file.exists('CODAindex.txt')){
		# This skips all of the output moving and checking - either the simulation has already been imported, or it was created using a custom method that put CODAindex.txt (and therefore hopefully everything else) in the root folder

		# end.time should have been re-saved along with the sim.chains by a previous readin, but just in case:
		if(!'end.time'%in%objects()){
			fi <- file.info(file.path(paste('sim',1:n.sims,sep='.'),'jagsoutput.txt'))
			end.time <- max(fi[,'mtime'])
		}
		
		visiblechains <- grep('CODAchain[[:digit:]+]', list.files(), value=TRUE)
		visnums <- as.numeric(gsub('.txt','',gsub('CODAchain','',visiblechains),fixed=TRUE))
		if(is.logical(sub.chains)){
			successful <- 1:n.chains %in% visnums
		}else{
			successful <- sub.chains %in% visnums
		}

        # Stop if all sims have crashed or we don't want a partial return:
		if(!all(successful)){
			pluralsim <- sum(!successful) > 1
			pluralchain <- length(unlist(sim.chains[!successful])) > 1
			if(n.sims==1)
				stop("The simulation appears to have crashed - check the model output in failed.jags() for clues", call.=FALSE)
			if(all(!successful))
				stop("All the simulations appear to have crashed - check the model output in failed.jags() for clues", call.=FALSE)
			if(!part.return)
				stop(paste("Simulation number", if(pluralsim) "s", " ", paste(which(!successful),collapse=", "), " appear", if(!pluralsim) "s", " to have crashed; check the output in failed.jags() for clues", sep=""), call.=FALSE)
		
			# Otherwise work out which chains correspond to these sims and just bin them unless we have been told to ignore them:
			if(any(unlist(sim.chains[!successful]) %in% which.chains))
				warning(paste("Chain", if(pluralchain) "s", " ", paste(unlist(sim.chains[!successful]),collapse=", "), " (simulation", if(pluralsim) "s", " ", paste(which(!successful),collapse=", "), ") appear", if(!pluralchain) "s", " to have crashed and will not be imported; check the output in failed.jags() for clues", sep=""))
			sim.chains <- sim.chains[successful]
			which.chains <- unlist(sim.chains)
			n.chains <- length(which.chains)
			simfolders <- paste("sim.",(1:n.sims)[successful],sep="")
			usingsims <- (1:n.sims)[successful]
			n.sims <- sum(successful)
		
			# Make sure the outputs are saved:
			forcefailed <- TRUE
		}
		
	}else{
	
		outputs <- vector('list', length=n.sims)
		for(i in 1:n.sims){
			tries <- 0
			repeat{
				outputs[[i]] <- try(paste(readLines(paste("sim.",i,"/jagsoutput.txt",sep=""), warn=FALSE), collapse="\n"), silent=TRUE)
				if(class(outputs[[i]])!="try-error")
					break

				#	Makes results.jags very slow:
#					simfinished <- grepl("[[:space:]]Deleting model",outputs[[i]])
					# Wait up to 3 seconds for 'Deleting model' to be appended:
#					if(simfinished || tries > 6)
#						break
#					if(runjags.getOption('debug'))
#						swcat('Failed to find "Deleting model" in model output - re-trying...\n')
#				}
				
				tries <- tries +1
				if(tries==11){
					if(n.sims==1) stop("An unexpected error occured: unable to read the output of the simulation") else stop(paste("An unexpected error occured: unable to read the output of simulation ", i, sep=""))
				}
				Sys.sleep(0.5)
			}
		}
		
		# If the model did not adapt, it will either have stopped and said Adaptation incomplete (adapt 0) or kept going and given a warning (update 0):
		notadapted <- any(sapply(outputs,function(x) return(grepl("[[:space:]]Adaptation incomplete",x))))
		# If the model required no adapt, it will either say Error (soon to be changed to note) .... not in adaptive mode  - or update 0 is ignored
		adaptunnecessary <- any(sapply(outputs,function(x) return(grepl("[[:space:]]not in adaptive mode",x))))
		stopifnot(!all(adaptunnecessary, notadapted))  # It should not be possible to get both messages
		# Check adaptation was completed
	
		if(target.adapt>0){
			# After the model is compiled, adapt(0) will either return an error or TRUE
			if(notadapted){	
				if(runjags.getOption('adapt.incomplete')=='error'){
					# Need to somehow kill rogue processes here but the PID information is lost by this point ...?					
					return(list(errormessage=paste('The adaptation phase of ', if(n.sims>1) 'one or more models ' else 'the model ', 'was not completed in ', target.adapt, ' iterations - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE))
				}
				if(runjags.getOption('adapt.incomplete')=='warning')
					warning(paste('The adaptation phase of ', if(n.sims>1) 'one or more models ' else 'the model ', 'was not completed in ', target.adapt, ' iterations, so the current samples may not be optimal - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
			}
			if(adaptunnecessary){
				if(!silent.jags && runjags.getOption('adapt.incomplete')!='silent')
					swcat('Note: the model did not require adaptation\n')
			}		
		}else{
			if(notadapted){
				if(runjags.getOption('adapt.incomplete')=='error')
					return(list(errormessage=paste('The ', if(n.sims>1) 'models require' else 'model requires', ' an adaptation phase (or longer burnin phase) - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE))
				if(runjags.getOption('adapt.incomplete')=='warning')
					warning(paste('The ', if(n.sims>1) 'models require' else 'model requires', ' an adaptation phase (or longer burnin phase), so the current samples may not be optimal - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
			}
		}	

		syntaxerrors <- sapply(outputs,function(x) return(grepl("[[:space:]]syntax error",x)))
		if(any(syntaxerrors)){
			scripts <- vector('list', length=n.sims)
			for(i in 1:n.sims){
				scripts[[i]] <- try(paste(readLines(paste("sim.",i,"/script.cmd",sep=""), warn=FALSE), collapse="\n"), silent=TRUE)
			}
			cat('Error: JAGS reported syntax errors in the command script: please send the following error message to the package mainainter:\n\n', paste(paste(outputs,scripts,sep="\n   -->\n"),collapse='\n\n'),'\n', sep='')
			stop('Error: JAGS reported syntax errors in the command script: please send the above error message to the package mainainter', call.=FALSE)
		}
		
		# Check that all simulations have finished (again - just to be sure):
		finished <- sapply(outputs,function(x) return(grepl("[[:space:]]Deleting model",x)))
		if(!all(finished)){
			if(suspended){
				allok <- TRUE
				unfinished <- TRUE
				return(list(unfinished=TRUE, finished.chains=finished))
			}else{
				stop("An unknown error occured - the simulation(s) appear to have not finished; check the model output in failed.jags() for clues")
			}
		}
			
		# Check for warnings about no monitored nodes ## MODIFY WHEN ALLOWING NO MONITORS ##:
		nomons <- sapply(outputs,function(x) return(grepl("[[:space:]]There are no monitors[[:space:]]",x)))
		if(any(nomons)) return(list(errormessage="The monitored nodes indicated do not exist in the model"))
		
		# Check that all simulations didn't crash - if batch mode there will only be an Updating 0 if successful, if not batch mode there will always be a can't update no model if unsuccessful:
		successful <- sapply(outputs,function(x){
			return(grepl("[[:space:]]Updating 0[[:space:]]",x) & !grepl("[[:space:]]Can't update. No model![[:space:]]",x))
			})
	
		usingsims <- (1:n.sims)
		if(!all(successful)){
			#### LOOK FOR THE CRASHED DUMP FILES HERE??? ####

			# Stop if all sims have crashed or we don't want a partial return:
			pluralsim <- sum(!successful) > 1
			pluralchain <- length(unlist(sim.chains[!successful])) > 1
			if(n.sims==1)
				stop("The simulation appears to have crashed - check the model output in failed.jags() for clues")
			if(all(!successful))
				stop("All the simulations appear to have crashed - check the model output in failed.jags() for clues")
			if(!part.return)
				stop(paste("Simulation number", if(pluralsim) "s", " ", paste(which(!successful),collapse=", "), " appear", if(!pluralsim) "s", " to have crashed; check the output in failed.jags() for clues", sep=""))
		
			# Otherwise work out which chains correspond to these sims and just bin them unless we have been told to ignore them:
			if(any(unlist(sim.chains[!successful]) %in% which.chains))
				warning(paste("Chain", if(pluralchain) "s", " ", paste(unlist(sim.chains[!successful]),collapse=", "), " (simulation", if(pluralsim) "s", " ", paste(which(!successful),collapse=", "), ") appear", if(!pluralchain) "s", " to have crashed and will not be imported; check the output in failed.jags() for clues", sep=""))
			
			sim.chains <- sim.chains[successful]
			which.chains <- unlist(sim.chains)
			n.chains <- length(which.chains)
			simfolders <- paste("sim.",(1:n.sims)[successful],sep="")
			usingsims <- (1:n.sims)[successful]
			n.sims <- sum(successful)
			
			# Make sure the outputs are saved:
			forcefailed <- TRUE
		} 
		
		# Now we have established all the simulations exited successfully, but there may be delays in obtaining the coda output...		
		
		# Modification dates may be wrong by max 2 seconds on Windows, hopefully less on unix
		mod.delay <- if(.Platform$OS.type=='unix') 1.5 else 3.5
		
		increment <- 0.5
		msgtime <- (mod.delay+2)/increment
		timeout <- runjags.getOption('timeout.import')/increment
		filesizeprop <- if(suspended) 0.1 else 0.5
		
		# First wait for all of the codaindex files to appear - this should be very fast:
		tries <- 0
		repeat{
			indexok <- file.exists(paste(simfolders,.Platform$file.sep,"CODAindex.txt",sep=""))
			if(all(indexok)) break
			tries <- tries +1
			if(tries==msgtime) swcat("Waiting for the CODA index files to be created...\n")
			if(tries==timeout) stop(paste("An unexpected error occured:  timed out waiting for the CODA index files to be created. The files available at time out were: ", paste(list.files(recursive=TRUE),collapse=", "), ". Please file a bug report (including this message) to the runjags package author.", sep=""))
			Sys.sleep(increment)
		}
		
		# Now wait for all of the codachain files to appear:
		codapaths <- unlist(lapply(1:n.sims,function(x) return(paste(simfolders[x],.Platform$file.sep,"CODAchain",1:length(sim.chains[[x]]),".txt",sep=""))))
		tries <- 0
		repeat{
			indexok <- file.exists(codapaths)
			if(all(indexok)) break
			tries <- tries +1
			if(tries==msgtime) swcat("Waiting for the CODA files to be created...\n")
			if(tries==timeout) stop(paste("An unexpected error occured:  timed out waiting for the CODA files to be created. The files available at time out were: ", paste(list.files(recursive=TRUE),collapse=", "), ". Please file a bug report (including this message) to the runjags package author.", sep=""))
			Sys.sleep(increment)
		}
		
		# Now wait for all of the codachain files to be finished writing out
		tries <- 0
		warninggiven <- FALSE
		repeat{

			fi <- file.info(codapaths)

			# Run a check on absolute file size here - give a note if > 1GB for any chain:
			if(!warninggiven && any(fi[,'size']>10^9)){
				swcat('NOTE: The JAGS output file(s) appear(s) to be very large - they may take some time to read.  Have you accidentally included a large vector in "monitor", or are you trying to run too many iterations without specifying "thin"?')
				swcat('If the read-in process fails (or is aborted), use ?results.jags and the read.monitor argument to retrieve the simulation')
				warninggiven <- TRUE
			}
 
			if(runjags.getOption('debug')>=5){
				swcat('Last modified (seconds ago) / file size as proportion of largest:\n')
				print(as.numeric(difftime(Sys.time(), fi[,'mtime'], units='secs')))
				print(fi[,'size']/max(fi[,'size']))
			}
			
			# Checking based purely on file sizes is not very sensitive but pretty specific and very fast when it works:
			if(n.chains>1 && all(fi[,'size']>0) && all(fi[,'size']/max(fi[,'size']) > 0.98)) break

			# Checking based on modification times with a small margin for file size is better (and the only possible method for 1 chain)
			if(all(fi[,'size']>0) && all(fi[,'size']/max(fi[,'size']) > filesizeprop) && all((fi[,'mtime'] + mod.delay) < Sys.time())) break
			
			tries <- tries +1
			if(tries==msgtime && !warninggiven) swcat("Waiting for the CODA files to be completed...\n")
				# If the coda files exist, wait a LONG time for them to be the correct sizes before giving up!
			if(tries==timeout) stop(paste("Timed out waiting for the CODA files to be completed. You can wait for the files to be written, and then run results.jags(\"", temp.directory, "\")\nThe file size and modification times at ", Sys.time(), " were: ", paste(paste(codapaths,fi[,'size'],as.character(fi[,'mtime']), sep=" : "),collapse=", "), ".  Please file a bug report (including this message) to the runjags package author.", sep=""))
			Sys.sleep(increment)
		}		
		end.time <- max(fi[,'mtime'])
		
		# At this point everything should be OK to proceed....
		
		for(s in 1:length(simfolders)){

			setwd(simfolders[s])
			
			# Try to copy the script file back to the main folder, mostly for my benefit
			# NO LONGER REQUIRED - we keep all folders
			# suppressWarnings(try(file.rename("script.cmd", paste("../script.",usingsims[s],".cmd",sep="")), silent=TRUE))
		
			# Copy the chains to the root simulation directory and renumber at the same time:
			fakenums <- 1:length(sim.chains[[s]])
			realnums <- sim.chains[[s]]
			for(c in fakenums){
				
				# One last repeat failsafe...
				tries <- 0
				repeat{
					success <- file.rename(paste("CODAchain",c,".txt",sep=""), paste("../CODAchain",realnums[c],".txt",sep=""))
					if(success) break
					tries <- tries +1
					if(tries==5) stop(paste("There was an error moving the coda file for chain ", realnums[c], sep=""))
					Sys.sleep(1)
				}
				close(file(paste("../CODAchain",realnums[c],".txt",sep="")))
			}

			# These were created before the coda files, so just assume they are OK:
			fm <- list.files(pattern="pd[[:graph:]]*.txt")
			file.rename(from=fm, to=file.path('..',fm))
			for(f in file.path("../",fm)) close(file(f))
			fm <- list.files(pattern="popt[[:graph:]]*.txt")
			file.rename(from=fm, to=file.path('..',fm))
			fm <- list.files(pattern="deviance[[:graph:]]*.txt")
			file.rename(from=fm, to=file.path('..',fm))
			for(f in file.path("../",fm)) close(file(f))

			if(s==1){
				success <- file.rename(from='CODAindex.txt', to='../CODAindex.txt')
				if(!success) stop("There was an error copying the coda index file")
				close(file("../CODAindex.txt"))
				for(f in c("CODAchain0.txt","CODAindex0.txt","CODAtable0.txt")){
					if(file.exists(f)){
						file.rename(from=f, to=file.path("..",f))
						close(file(file.path("..",f)))
					}
				}
				suppressWarnings(success <- try({
					fs <- read.csv('samplers.csv',sep='\t',header=FALSE,stringsAsFactors=FALSE)
					success <- file.rename(from='samplers.csv', to='../samplers.csv')
					close(file("../samplers.csv"))
        		}, silent=TRUE))
        		if(inherits(success, 'try-error'))
          		  fs <- NA
			}
			
			if(runjags.getOption('debug') && !identical(fs,NA) && s!=1){
			  suppressWarnings(try({
				ns <- (read.csv('samplers.csv',sep='\t',header=FALSE,stringsAsFactors=FALSE))
				if(runjags.getOption('debug') && (identical(ns,NA) || !all(ns==fs))){
					cat('SAMPLERS NOT EQUAL\n')
					browser()
				}
				}, silent=TRUE))
			}
      
			setwd(temp.directory)
			# Never delete the sim folder - it has command and potentially interesting stuff in it!!  Also we need to read the jagsoutput for catting to screen.
			# unlink(simfolders[s], recursive=TRUE)
		}
	}
	
	setwd(temp.directory)
	
	chainscopied <- sort(as.numeric(gsub("[CODAchain,.txt]", "", list.files(pattern="CODAchain[[:digit:]]+[.]txt"))))
	# We might have a codachain0 if monitoring DIC so remove from the fakenums:
	chainscopied <- chainscopied[chainscopied!=0]					

	if(length(chainscopied) < n.chains || (length(chainscopied) > n.chains) && !part.return){
		stop(paste("Expected ", n.chains, " chains to be output but found ", length(chainscopied), " in the root simulation directory - please file a bug report to the runjags package author", sep=""))		
	}
	
	allok <- TRUE


	swcat("Simulation complete.  Reading coda files...\n")

	
	# Make sure sub samples makes sense - convert to NA or a number:
	if(identical(sub.samples, FALSE)){
		sub.samples <- NA
	}else{
		if(sub.samples > target.iters){
			warning('Specified sub.samples was greater than the number of iterations and was ignored', call.=FALSE)
			sub.samples <- NA
		}else{
			target.iters <- sub.samples
		}
	}
	
	suppressWarnings(inputsuccess <- try(input.data <- read.openbugs.subset(start=NA, end=NA, thin=NA, quiet=TRUE, vars=read.monitor, sub.samples=sub.samples, which.chains=which.chains), silent=TRUE))
	if((class(inputsuccess)=="try-error")){
		filename <- paste("jags.dump", 1, ".R", sep="")
		suppressWarnings(try(readlinessuccess <- try(tempinput <- readLines(filename)), silent=TRUE))
		if(class(readlinessuccess)=="try-error"){
			
			stop(paste("Unable to load model output - the following error was returned from read.openbugs.subset:  ", as.character(inputsuccess), ". Please file a bug report to the runjags package author.", sep=""))
		}else{
			stop("The model appears to have crashed during the burnin period.  Check failed.jags() for clues.")
		}
	}
	
	# Last failsafe:
	if(nchain(input.data)!=n.chains) stop(paste("Expected ", n.chains, " chains to be returned but only found ", nchain(input.data), " - please file a bug report to the runjags package author", sep=""))
	
	swcat("Coda files loaded successfully\n")
	achieved <- niter(input.data)

	if(any(c("full.pd") %in% monitor)){
		pd <- try(read.coda.subset('CODAchain0.txt','CODAindex0.txt', start=NA, end=NA, thin=NA, quiet=TRUE, vars=read.monitor, sub.samples=sub.samples), silent=TRUE)
		if(class(pd)=="try-error"){
			swcat("\n")			
			warning("There was an error reading the full sum(pD)")
		}
	}else{
		pd <- NA
	}
	if(any(monitor=="dic")){
		dn <- ''
		pdtab <- try(read.table('pdtable0.txt', header=FALSE), silent=TRUE)
		if(class(pdtab)=="try-error"){
			warning("There was an error reading the mean pD table")
			mean.pd <- NA
		}else{
			mean.pd <- as.matrix(pdtab[,2])
			dn <- pdtab[,1]
		}
		popttab <- try(read.table('popttable0.txt', header=FALSE), silent=TRUE)
		if(class(popttab)=="try-error"){
			warning("There was an error reading the mean pOpt table")
			mean.popt <- NA
		}else{
			mean.popt <- as.matrix(popttab[,2])
			dn <- popttab[,1]
		}
		devtab <- try(read.table('deviancetable0.txt', header=FALSE), silent=TRUE)
		if(class(devtab)=="try-error"){
			warning("There was an error reading the mean deviance table")
			mean.deviance <- NA
		}else{
			mean.deviance <- as.matrix(devtab[,2])
			dn <- devtab[,1]
		}		
		
		deviance.table <- cbind(mean.deviance, mean.pd, mean.popt)
		dimnames(deviance.table) <- list(dn, c('mean.deviance', 'mean.pD', 'mean.pOpt'))
		deviance.sum <- apply(deviance.table,2,sum)
		names(deviance.sum) <- c('sum.mean.deviance', 'sum.mean.pD', 'sum.mean.pOpt')
	}else{
		deviance.table <- NA
		deviance.sum <- NA
	}
	
	suppressWarnings(success <- try({
		samplers <- read.csv('samplers.csv',sep='\t',header=FALSE,stringsAsFactors=FALSE)
		dimnames(samplers)[[2]] <- c('Index.No', 'Sampler', 'Node')
		samplers$Index.No <- as.numeric(samplers$Index.No)
	}, silent=TRUE))
	if(inherits(success, 'try-error'))
	  samplers <- NA

	if(!is.na(sub.samples) && abs(achieved-target.iters)>0.1){
		crashed <- TRUE
		swcat("Warning:  Simulation crashed after ", achieved, " iterations (expecting ", target.iters, ")\n", sep="")
		
		crash.end <- character(length=n.chains)
		for(i in 1:n.chains){
			filename <- paste("jags.dump", i, ".R", sep="")
			suppressWarnings(inputsuccess <- try(tempinput <- readLines(filename), silent=TRUE))
			if(class(inputsuccess)=="try-error"){
				swcat("Error reading crash point of chain ", i, "\n", sep="")
				crash.end[i] <- NA
			}else{
				crash.end[i] <- ""
				for(j in 1:length(tempinput)){
					crash.end[i] <- paste(crash.end[i], tempinput[j], "\n", sep="")
				}
			}
		}
	}else{
		crashed <- FALSE
		input.end <- character(length=n.chains)
		for(i in 1:n.chains){
			filename <- paste("out", which.chains[i], ".Rdump", sep="")
			suppressWarnings(inputsuccess <- try(tempinput <- readLines(filename), silent=TRUE))
			if(class(inputsuccess)=="try-error"){
				swcat("Error reading end point of chain ", which.chains[i], ".\n", sep="")
				input.end[i] <- NA
			}else{
				input.end[i] <- ""
				for(j in 1:length(tempinput)){
					input.end[i] <- paste(input.end[i], tempinput[j], "\n", sep="")
				}
			}
		}
	}
	



	if(any(is.na(unlist(input.data)))){
	
		nastring <- ""
	
		varnames <- dimnames(input.data[[1]])[[2]]
	
		for(i in 1:n.chains){
		
			for(j in 1:nvar(input.data)){
			
				if(any(is.na(input.data[[i]][,j]))){
				
					nastring <- paste(nastring, if(nastring!="") ", ", varnames[j], " (chain ", i, ")", sep="")
				
				}
			
			}
		
		}
		stop(paste("One or more of the values for the monitored variable(s) '", nastring, "' was invalid (missing data).  Ensure that the model syntax is correct and that appropriate starting values have been given.", sep=""))
	}

	if(crashed) inits <- unlist(crash.end) else inits <- unlist(input.end)
	class(inits) <- 'runjagsinits'
	
	# For future re-imports:
	sim.chains <- orig.sim.chains
	save(sim.chains, end.time, file='simchainsinfo.Rsave')
	
	return(list(mcmc=input.data, deviance.table=deviance.table, deviance.sum=deviance.sum, pd=pd, end.state=inits, samplers=samplers, end.time=end.time))
}



############################
# The following two functions were modified (in October 2014) from functions in coda version 0.16-1
# Modifications add the possibility of being selective about which variables to import (in case of large coda files)
# Some other modifications regarding missing->NA variables were also required to work with runjags, and the bug preventing 1 iteration being imported was fixed
# Original functions were GPL>=2, copyright of Martyn Plummer
# These functions are not exposed and will be removed when an updated version of coda with similar modifications becomes available
# Note: checkvalidmonitorname and expandindexnames are currently required for read.coda.subset
############################

read.openbugs.subset <- function (stem = "", start, end, thin, quiet = FALSE, vars, sub.samples, which.chains)
	{
		# runjags will pass NA not missing:
		if(missing(start))
			start <- NA
		if(missing(end))
			end <- NA
		if(missing(thin))
			thin <- NA
		if(missing(vars))
			vars <- NA
		if(missing(sub.samples))
			sub.samples <- NA
		if(missing(which.chains))
			which.chains <- NA
			
	    index.file <- paste(stem, "CODAindex.txt", sep = "")
	    if (!file.exists(index.file)) 
	        stop("No index file found")
	    index.date <- file.info(index.file)$ctime

		# Sub select chains:
		if(!identical(which.chains, NA) && length(which.chains)>0){
	        output.file <- paste(stem, "CODAchain", which.chains, ".txt", sep = "")
			if(any(!file.exists(output.file))){
				warning('The which.chains argument specified more than the available number of chains', call.=FALSE)
				which.chains <- which.chains[file.exists(output.file)]
			}
			output.file <- output.file[which.chains]
			nchain <- length(which.chains)
		}else{
		    nchain <- 0
		    while (TRUE) {
		        output.file <- paste(stem, "CODAchain", nchain + 1, ".txt", 
		            sep = "")
		        if (file.exists(output.file)) {
		            nchain <- nchain + 1
		            output.date <- file.info(output.file)$ctime
		            dt <- difftime(index.date, output.date, units = "mins")
		            if (abs(as.numeric(dt)) > 1) {
		                warning(paste("Files \"", index.file, "\" and \"", 
		                  output.file, "\" were created at different times\n", 
		                  sep = ""))
		            }
		        }
		        else break
		    }
			which.chains <- 1:nchain
		}
		
	    if (nchain == 0) 
	        stop("No output files found")
		
	    ans <- vector("list", nchain)
	    for (i in 1:nchain) {
	        output.file <- paste(stem, "CODAchain", which.chains[i], ".txt", sep = "")
	        ans[[i]] <- read.coda.subset(output.file, index.file, start, 
	            end, thin, quiet, vars, sub.samples)
	    }

		return(mcmc.list(ans))
	}



read.coda.subset <- function (output.file, index.file, start, end, thin, quiet = FALSE, vars, sub.samples) 
{
	# runjags will pass NA not missing:
	if(missing(start))
		start <- NA
	if(missing(end))
		end <- NA
	if(missing(thin))
		thin <- NA
	if(missing(sub.samples))
		sub.samples <- NA
	
	if(!is.na(sub.samples) && sub.samples < 1) stop('sub.samples must be >=1')
	if(missing(vars)) vars <- NA else vars <- checkvalidmonitorname(vars)
	
	index <- read.table(index.file, row.names = 1, col.names = c("", "begin", "end"))
    vnames <- row.names(index)
	
	iterations <- index[1,2]
	# TODO:  Modify later
	
    if (is.R()) {
      	
		if(!identical(as.character(vars), as.character(NA))){
			
			# Indexes of vars will be expanded earlier in runjags code
			# For variables provided non-indexed:
			nonindvars <- paste(vars,'\\[.*\\]',sep='')

			matchvarsfun <- function(x) return(any(x==vars) || any(sapply(nonindvars,grepl,x=x)))
			readin <- which(sapply(vnames,matchvarsfun))
			if(length(readin)==0)
				stop(paste('None of the specified monitor(s) "', paste(vars, collapse='", "'), '" were found in the coda files', sep=''), call.=FALSE)
			vars <- vnames[readin]
		
			starts <- index[readin[1],1]
			ends <- numeric(0)
			if(length(readin) > 1){
				for(i in 2:length(readin)){
				  if(index[readin[i],1]!=(index[readin[i-1],2]+1)){
				    ends <- c(ends, index[readin[i-1],2])
				    starts <- c(starts, index[readin[i],1])  
				  }
				  last <- index[readin[i],2]
				}
			}
			ends <- c(ends, index[readin[length(readin)],2])

			vnames <- vnames[readin]
			
			temp <- as.data.frame(matrix(as.numeric(NA), nrow=iterations*length(vnames), ncol=2, dimnames=list(NULL, c('iter', 'val'))))
			
			ind.s <- 1
			for(f in 1:length(starts)){
				ind.e <- ind.s + (ends[f]-starts[f])
				temp[ind.s:ind.e,] <- as.data.frame(scan(output.file, what = list(iter = 0, val = 0), nmax=((ends[f]-starts[f])+1), skip=(starts[f]-1), quiet = TRUE))
        ind.s <- ind.e+1
				# TODO: index scan to immediately drop thinned chains ... make sure pd thinned later
			}
			
			index <- index[readin,,drop=FALSE]
			index[,1] <- (0:(length(vnames)-1) * iterations)+1
			index[,2] <- 1:length(vnames) * iterations
			
		}else{
   		 	# as.data.frame prevents 1 iteration and 1 variable being a list
	        temp <- as.data.frame(scan(output.file, what = list(iter = 0, val = 0), quiet = TRUE))
			# TODO: index scan to immediately drop thinned chains ... make sure pd thinned later
		}

		
    }
    else {
		stop('Modified function is not S-compatible')
		# NOT modified as I don't use S...
		if(!identical(vars, NA)) warning("Argument 'vars' was ignored")
        temp <- scan(output.file, what = list(iter = 0, val = 0))
		# TODO: index scan to immediately drop thinned chains ... make sure pd thinned later
    }
	
	# Added to allow sub.samples to control thinning indirectly:
	if(!is.na(sub.samples)){
		if(sub.samples > iterations) warning('Specified sub.samples was greater than the number of iterations and was ignored', call.=FALSE)
		thin <- max(1,floor(iterations/sub.samples))
	}
	
	  start.vec <- end.vec <- thin.vec <- numeric(nrow(index))
    for (i in 1:length(vnames)) {
        iter.i <- temp$iter[index[i, "begin"]:index[i, "end"]]
		
		# Modified to allow a single iteration to be read:
        thin.i <- max(1, unique(diff(iter.i)))
        thin.vec[i] <- if (length(thin.i) == 1) thin.i else NA
        start.vec[i] <- iter.i[1]
        end.vec[i] <- iter.i[length(iter.i)]
    }
    if (any(is.na(start.vec)) || any(thin.vec != thin.vec[1]) || 
        any((start.vec - start.vec[1])%%thin.vec[1] != 0)) {
        iter <- sort(unique(temp$iter))
        old.thin <- unique(diff(iter))
        if (length(old.thin) == 1) 
            is.regular <- TRUE
        else {
            if (all(old.thin%%min(old.thin) == 0)) 
                old.thin <- min(old.thin)
            else old.thin <- 1
            is.regular <- FALSE
        }
    }
    else {
        iter <- seq(from = min(start.vec), to = max(end.vec), 
            by = thin.vec[1])
        old.thin <- thin.vec[1]
        is.regular <- TRUE
    }
    if (is.na(start)) 
        start <- min(start.vec)
    else if (start < min(start.vec)) {
        warning("start not changed")
        start <- min(start.vec)
    }
    else if (start > max(end.vec)) 
        stop("Start after end of data")
    else iter <- iter[iter >= start]
    if (is.na(end)) 
        end <- max(end.vec)
    else if (end > max(end.vec)) {
        warning("end not changed")
        end <- max(end.vec)
    }
    else if (end < min(start.vec)) 
        stop("End before start of data")
    else iter <- iter[iter <= end]
		
    if (is.na(thin)) 
        thin <- old.thin
    else if (thin%%old.thin != 0) {
        thin <- old.thin
        warning("thin not changed")
    }
    else {
        new.iter <- iter[(iter - start)%%thin == 0]
        new.thin <- unique(diff(new.iter))
        if (length(new.thin) != 1 || new.thin != thin) 
            warning("thin not changed")
        else {
            iter <- new.iter
            end <- max(iter)
            is.regular <- TRUE
        }
    }
    out <- matrix(NA, nrow = length(iter), ncol = nrow(index))
    dimnames(out) <- list(iter, vnames)
    for (v in vnames) {
        if (!quiet) 
            swcat("Abstracting", v, "... ")
        inset <- index[v, "begin"]:index[v, "end"]
        iter.v <- temp$iter[inset]
        if (!is.regular) {
            use.v <- duplicated(c(iter, iter.v))[-(1:length(iter))]
            use <- duplicated(c(iter.v, iter))[-(1:length(iter.v))]
        }
        else {
            use.v <- (iter.v - start)%%thin == 0 & iter.v >= 
                start & iter.v <= end
            use <- (iter.v[use.v] - start)%/%thin + 1
        }
        if (length(use) > 0 && any(use.v)) 
            out[use, v] <- temp$val[inset[use.v]]
        if (!quiet) 
            swcat(length(use), "valid values\n")
    }
    if (is.regular) 
        out <- mcmc(out, start = start, end = end, thin = thin)
    else warning("not returning an mcmc object")
		
	# Added to allow sub.samples to control thinning indirectly:
	if(!is.na(sub.samples)){
		out <- as.mcmc(out[1:sub.samples,,drop=FALSE])
	}	
	
    return(out)
}

############################
