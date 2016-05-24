# Methods can return either a character string representing an error to be returned, or TRUE if files are to be read immediately or FALSE if JAGS run is ongoing.  Or a list starting with TRUE or FALSE (optionally named complete), and containing other things that will be returned in the runjags object if the JAGS run is ongoing.  If it returns an error the JAGS run is assumed to be not successful.

runjags.simple <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, rjags){
	
	retval <- 'An unknown error occured while calling JAGS using the simple method'
	
	stopifnot(n.sims==1)
	
	if(runjags.getOption('debug'))
		swcat('Using the following simulation directory: ', getwd(), '\n')
	
	tryCatch({

		swcat("Running the simulation using the simple method... (output will be displayed once the simulation has termianted)\n")
		flush.console()
  	  	
			# In theory more portable code although can't guarantee interleaving of stdout/stderr which BREAKS run.jags.study:
			# success <- system2(jags, args=if(batch.jags) "sim.1/script.cmd" else character(0), stdout="", stderr="", stdin=if(batch.jags) "" else "sim.1/script.cmd")
			
    	if(os == "windows"){
			# This was working but sometimes seems not to fully quote the path to JAGS?  One above doesn't work either...
			suppressWarnings(success <- shell(paste(shQuote(jags, type='cmd'), if(!batch.jags) " <", " sim.1/script.cmd 2>&1", sep = ""), intern=TRUE, wait=TRUE, ignore.stderr=FALSE))

		}else{
			suppressWarnings(success <- system(paste(shQuote(jags), if(!batch.jags) " <", " sim.1/script.cmd 2>&1", sep=""), intern=TRUE, wait=TRUE, ignore.stderr=FALSE))
		}
	

		if(!silent.jags){
			cat(success,sep='\n')
		}
		cat(success,file="sim.1/jagsoutput.txt",sep="\n")
		
		flush.console()
		
		retval <- TRUE
		
	})

	return(retval)
}

runjags.rjparallel <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, rjags){

	# Eventually we should compile the model once and then export the compiled model and then change the inits - but we are currently unable to change RNG states on compiled models with rjags 3.x
	# We also have to destroy the compiled rjags object after calling runjags.start() as we update different compiled objects

	# rjags object should have been compiled by extend.jags or autoextend.jags if this was necessary - and the modules loaded just before that
	# BUT will have to be re-compiled with appropriate RNG on each node
	# but having already compiled it we have at least multiple safe RNG streams set up that we can just copy...

	extra.options <- rjags[names(rjags)!='rjags']
	rjags <- rjags$rjags
	
	# There should never be a pd, popt or full.pd monitor here but remove just in case:
	monitor <- extra.options$monitor[! extra.options$monitor %in% c("pd","popt","full.pd","dic","ped")]
	
	# Should be ignoring by and progress.bar - always none
#	by <- if(is.na(extra.options$by)) min(100, extra.options$burnin/50) else extra.options$by
	
	# Retrieve model, data and initial values with RNG states (guaranteed to be thread-proof by auto/extend.jags code) for all chains:
	model <- paste(rjags$model(),collapse="\n")
	data <- rjags$data()
	
	# Prevents problems with data/inits reloading:
	data <- data[names(data) %in% extra.options$origdatanames]
	
	inits <- rjags$state(internal=TRUE)
	# This should have .RNG.name and .RNG.state always:
	if(!all(sapply(inits,function(x) all(c('.RNG.name','.RNG.state')%in%names(x))))){
		if(runjags.getOption('debug')){
			cat('Inits .RNG.name extraction error...\n')
			browser()
		}
		stop('An error occured while re-extracting the initial values from the rjags object - you should be able to avoid this by specifying .RNG.name in the initial values', call.=FALSE)
	}		
	
	# If we have loaded a module externally then it needs to be loaded on the cluster too
	# The obvious one is lecuyer - if the model was compiled with it, we need to add it:
	if(any(grepl("lecuyer::RngStream",inits))){
		extra.options$modules <- checkmodfact(c(extra.options$modules, list(c("lecuyer","TRUE"))), 'module')
	}	
	
	outfile <- ''
	if(identical(cl, NA)){
		outfile <- tempfile()
		cl <- parallel::makeCluster(n.sims, type=if(.Platform$OS.type=='unix') 'FORK' else 'PSOCK', outfile=outfile)
		on.exit({
			stopCluster(cl)
			unlink(outfile)
		})
	}
	
	clname <- paste(summary(cl)[1,'Class'], capture.output(print(cl)))
	clname <- gsub('forknode socket', 'Fork', clname)
	clname <- gsub('SOCKnode socket', 'PSOCK', clname)
	
	# We can check if the model is adapted here (don't actually do anything, just check):
	if(checkadaptrequired(rjags)){
		if(extra.options$adapt==0 && extra.options$burnin <1000){   # If burnin is long enough, assume that adaptation is OK until I can do adapt(,0)
			if(runjags.getOption('adapt.incomplete')=='error')
				stop(paste('The model requires an adaptation phase and the burnin period is not long enough to guarantee adaptation - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
			if(runjags.getOption('adapt.incomplete')=='warning')
				warning(paste('The model requires an adaptation phase and the burnin period is not long enough to guarantee adaptation, so the current samples may not be optimal - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
		}	
	}else{
		if(extra.options$adapt>0 && !silent.jags && runjags.getOption('adapt.incomplete')!='silent')
			swcat('Note: the model did not require adaptation\n')
	}		
	
	# Only relevant if it is stop:
	extra.options$adapt.incomplete <- runjags.getOption('adapt.incomplete')
	
	if(n.sims==1) swcat("Starting the rjags simulation using a ", clname, "\n") else swcat("Starting ", n.sims, " rjags simulations using a ", clname, "\n")		
	
	clfun <- function(s,model,data,inits,extra.options,monitor){
		
		if(!require("runjags") || packageVersion('runjags')$major<2) stop("The runjags package (version >= 2) is not installed on one or more cluster nodes")
		if(!requireNamespace("rjags")) stop("The runjags package is not installed on one or more cluster nodes")
	
		# Module loading MUST be done before model compilation:
		if(!identical(extra.options$modules,"")) for(i in 1:length(extra.options$modules)){
			if(extra.options$modules[[i]][1]=="runjags"){
				if(extra.options$modules[[i]][2]=='TRUE'){
					success <- load.runjagsmodule()
				}else{
					success <- unload.runjagsmodule()
				}
			}else{
				if(extra.options$modules[[i]][2]=='TRUE'){
					success <- try(rjags::load.module(extra.options$modules[[i]][1]))
					if(is.null(success)) success <- TRUE
				}else{
					suppressWarnings(success <- try(rjags::unload.module(extra.options$modules[[i]][1])))
					if(is.null(success)) success <- TRUE
				}
			}
			if(inherits(success, 'try-error')) stop(paste("Failed to ", if(extra.options$modules[[i]][2]=='FALSE') "un", "load the module '", extra.options$modules[[i]][1], "'", sep=""))			
		}		
		if(!identical(extra.options$factories,"")) for(i in 1:length(extra.options$factories)){
			fa <- ""
			try(fa <- as.character(rjags::list.factories(extra.options$factories[[i]][2])$factory))
			if(!extra.options$factories[[i]][1] %in% fa) stop(paste("The factory '", extra.options$factories[[i]][1], "' of type '", extra.options$factories[[i]][2], "' is not available - ensure any required modules are also provided", sep=""))

			success <- try(rjags::set.factory(extra.options$factories[[i]][1], extra.options$factories[[i]][2], as.logical(extra.options$factories[[i]][3])))
			if(is.null(success)) success <- TRUE

			if(inherits(success, 'try-error')) stop(paste("Failed to ", if(extra.options$factories[[i]][3]=='FALSE') "un", "set the factory '", extra.options$factories[[i]][1], "' of type '", extra.options$factories[[i]][2], "'", sep=""))			
		}
		
		chains <- extra.options$sim.chains[[s]]
		inits <- inits[chains]
		n.chains <- length(chains)
		tmodel <- textConnection(model)
		
		rjags <- rjags::jags.model(file=tmodel, data=data, inits=inits, n.chains, n.adapt=0, quiet=TRUE)
		close(tmodel)
		
		adaptdone <- TRUE
		# If adapt=0 and burnin too small I have no way of telling if the model is adapted but rjags does print 'Stopping adaptation'... but ignore this and assume burnin >=1000 is OK
		flush.console()
		if(extra.options$adapt>0){
			adaptdone <- rjags::adapt(rjags,n.iter=extra.options$adapt,end.adaptation=TRUE,progress.bar='none')
			flush.console()
			if(!adaptdone && extra.options$adapt.incomplete=='error')
				return(list(adaptdone=adaptdone))
		}
		
		if(extra.options$burnin>0){
			flush.console()
			by <- if(is.na(extra.options$by)) min(100, extra.options$burnin/50) else extra.options$by
			update(rjags,n.iter=extra.options$burnin,progress.bar='none')
			flush.console()
		}	
		
		flush.console()
		samples <- rjags::jags.samples(rjags,variable.names=monitor,n.iter=extra.options$sample,progress.bar='none', thin=extra.options$thin)
		# Do this BEFORE our dummy call of 1 iteration
		end.state <- sapply(rjags$state(internal=TRUE),dump.format)	
			
		# This is just a dummy call so that we can get the names of the variables:
		suppressWarnings(varnames <- varnames(rjags::coda.samples(rjags,variable.names=monitor[monitor!="pD"],n.iter=2,progress.bar="none", thin=1)))
		flush.console()
	
		mcmcout <- lapply(samples, function(x){
			x <- coda::as.mcmc.list(x)
			names(x) <- as.character(chains)
			return(x)
		})
		
		nvar <- length(varnames)
		niter <- sapply(mcmcout,niter)
		if(!all(niter==niter[1])) stop("An error occured with the rjags method - variables returned with differing numbers of iterations")
		niter <- niter[1]
	
		mcmc <- vector('list',length=n.chains)
	
		# rjags (verison 3) returns variables in alphabetical order regardless of input order, so re-order here:
		reindex <- matchvars(monitor[monitor!="pD"],varnames,exactneeded=TRUE,testfound=FALSE)
		varnames <- varnames[reindex]
	
		for(i in 1:n.chains){
			mcmc[[i]] <- coda::mcmc(do.call('cbind',lapply(mcmcout, function(x) return(x[[i]])))[,reindex,drop=FALSE], start=extra.options$burnin+1, thin=extra.options$thin)
			dimnames(mcmc[[i]]) <- list(1:niter, varnames)
		}
		mcmc <- coda::as.mcmc.list(mcmc)
	
		
		success <- try({
			samplers <- rjags::list.samplers(rjags)
			samplers <- as.data.frame(t(matrix(unlist(lapply(1:length(samplers), function(x) return(rbind(Index.No=x, Sampler=names(samplers)[x], Node=samplers[[x]])))),nrow=3)), stringsAsFactors=FALSE)
			dimnames(samplers)[[2]] <- c('Index.No', 'Sampler', 'Node')
			samplers$Index.No <- as.numeric(samplers$Index.No)	
			}, silent=TRUE)
		if(inherits(success, 'try-error'))
			samplers <- NA
		
		return(list(mcmc=mcmc, end.state=end.state, samplers=samplers, adaptdone=adaptdone))
	
	}
	success <- try({
		allm <- parLapply(cl,1:n.sims,clfun,model=model,data=data,inits=inits,extra.options=extra.options,monitor=monitor)
		#	allm <- lapply(1:n.sims,clfun,model=model,data=data,inits=inits,extra.options=extra.options,monitor=monitor)
	}, silent=TRUE)
	if(inherits(success, 'try-error')){
		if(!file.exists(outfile))
			err <- paste("One or more rjags sessions failed with the following error:\n",as.character(success),'\nHave you remembered to specify all required modules and factories?',sep='')
		else
			err <- paste("One or more rjags sessions failed with the following error:\n",as.character(success),'\nThe worker log file (which may also help with debugging) is:\n',paste(readLines(outfile, warn=FALSE),collapse='\n'),'\nHave you remembered to specify all required modules and factories?',sep='')
		
		if(runjags.getOption('debug')){
			cat(err)
		}
		stop(err,call.=FALSE)
	}
	
	if(any(sapply(allm, function(x) return(!x$adaptdone)))){
		if(runjags.getOption('adapt.incomplete')=='error')
			stop(paste('The adaptation phase of one or more models was not completed in ', extra.options$adapt, ' iterations - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
		if(runjags.getOption('adapt.incomplete')=='warning')
			warning(paste('The adaptation phase of one or more models was not completed in ', extra.options$adapt, ' iterations, so the current samples may not be optimal - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
	}
	
	tvarnames <- sapply(allm,function(x) return(paste(varnames(x$mcmc),collapse=',')))
	if(!all(tvarnames==tvarnames[1])) stop("An error occured with the rjparallel method - simulations returned differing variable names")

	# Re-order and combine mcmc lists:	
	sim.chains <- extra.options$sim.chains
	mcmc <- as.mcmc.list(lapply(1:extra.options$n.chains, function(x){
		sim <- which(sapply(sim.chains,function(y) return(any(x==y))))
		return(allm[[sim]]$mcmc[[which(sim.chains[[sim]]==x)]])
	}))
	end.state <- sapply(1:extra.options$n.chains, function(x){
		sim <- which(sapply(sim.chains,function(y) return(any(x==y))))
		return(allm[[sim]]$end.state[[which(sim.chains[[sim]]==x)]])
	})
	
	samplers <- allm[[1]]$samplers
	if(runjags.getOption('debug') && !identical(NA, samplers)){
		allsamplers <- lapply(allm, function(x){
			t <- as.matrix(x$samplers)
			dimnames(t) <- list(x$samplers$Node,dimnames(t)[[2]])
			t <- t[samplers$Node,,drop=FALSE]
			return(t)
		})
    
		if(!all(sapply(allsamplers, function(x) return(all(as.numeric(x[,1])==samplers$Index.No) && all(x[,2]==samplers$Sampler))))){
			cat('SAMPLERS NOT EQUAL\n')
			browser()
		}		
	}
	
	if(!silent.jags) swcat("Simulation complete\n")
	
	return(list(complete=TRUE, mcmc=mcmc, deviance.table=NA, deviance.sum=NA, pd=NA, end.state=end.state, samplers=samplers))
	
}

runjags.snow <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, rjags){
		
	retval <- 'An unknown error occured while calling JAGS using the snow method'
	
	jags <- remote.jags
	
	files <- vector('list', length=n.sims+1)	
	initfiles <- list.files(pattern='inits')
	files[[1]] <- vector('list', length=length(initfiles)+2)
	names(files[[1]]) <- c(initfiles, 'model.txt', 'data.txt')
	for(i in 1:length(initfiles)){
		files[[1]][[i]] <- paste(readLines(initfiles[i], warn=FALSE), '\n', collapse="\n")
	}
	files[[1]]$model.txt <- paste(readLines("model.txt", warn=FALSE), '\n', collapse="\n")
	files[[1]]$data.txt <- paste(readLines("data.txt", warn=FALSE), '\n', collapse="\n")
	
	thed <- getwd()		
	for(s in 1:n.sims){
		
		setwd(paste("sim",s,sep="."))

		files[[s+1]] <- list()
		file <- list.files()
		for(f in 1:length(file)){
			files[[s+1]][file[f]] <- paste(readLines(file[f], warn=FALSE), '\n', collapse="\n")
		}
		
		setwd(thed)
	}
	
	tryCatch({

		f <- function(s, files, jags, batch.jags, silent.jags, path){
			
			if(!require("runjags")) return(paste("The runjags package was not found on the snow node '", Sys.info()['nodename'], "'", sep=""))

			retval <- "An error occured on the snow cluster"

			tryCatch({
							
				mdfiles <- files[[1]]
				simfiles <- files[[s+1]]
				if(s>1) silent.jags <- TRUE

				if(class(jags)=="function") jags <- jags()					
				if(jags=="*//*usefindjags*//*") jags <- findjags()

				testjags <- testjags(jags, silent=TRUE)
				if(!testjags$JAGS.available && !testjags$JAGS.found){
					return(paste("JAGS was not found on the snow node '", Sys.info()['nodename'], "' at the path '", jags, "'", sep=""))
				}
				jags <- testjags$JAGS.path
				
				setwd(path)
				existingfiles <- unique(c(list.files(), list.files(recursive=TRUE)))
				# Check to see that the sim files are there (which they will be if we created the cluster inside the function), if not create the folder:
				if(paste('sim',s,sep='.') %in% existingfiles){
					cleanup <- FALSE
				}else{
					cleanup <- TRUE
					
					success <- dir.create(paste('sim',s,sep='.'))
					newfiles <- unique(c(list.files(), list.files(recursive=TRUE)))

					if(!success)
						stop('An error occured while attempting to recreate the simulation folder')

					for(f in 1:length(mdfiles)){
						cat(mdfiles[[f]], '\n', file=names(mdfiles)[f])
					}
					setwd(paste("sim",s,sep="."))
					for(f in 1:length(simfiles)){
						cat(simfiles[[f]], '\n', file=names(simfiles)[f])
					}
				}
				
				setwd(path)

				newfiles <- unique(c(list.files(), list.files(recursive=TRUE)))
				if(any(!names(mdfiles) %in% newfiles))
					stop(paste('An error occured while attempting to recreate the model files\nPre-existing files:  ', paste(existingfiles,collapse=', '), '  -  after creating:  ', paste(newfiles,collapse=', '), sep=''))
				if(any(!file.exists(file.path(paste("sim",s,sep="."), names(simfiles)))))
					stop(paste('An error occured while attempting to recreate the simulation files\nPre-existing files:  ', paste(existingfiles,collapse=', '), '  -  after creating:  ', paste(newfiles,collapse=', '), sep=''))
				
				os <- .Platform$OS.type
			
				retval <- "An error occured while calling JAGS"
				
				if (os == "windows"){		
					currentsyspath <- Sys.getenv('PATH')
					if(!grepl(libpaths$PATH,currentsyspath,fixed=TRUE)){
						Sys.setenv(PATH=paste(currentsyspath, ';', testjags$libpaths$PATH, sep=''))
					}

					currentsysbinpath <- Sys.getenv('LTDL_LIBRARY_PATH')
					if(!grepl(libpaths$LTDL_LIBRARY_PATH,currentsysbinpath,fixed=TRUE)){
						Sys.setenv(LTDL_LIBRARY_PATH=paste(currentsysbinpath, if(currentsysbinpath!='') ';', testjags$libpaths$LTDL_LIBRARY_PATH, sep=''))
					}		
								
					success <- shell(paste(shQuote(jags), if(!batch.jags) " <", " sim.", s, "/script.cmd > sim.", s, "/jagsoutput.txt 2>&1", sep = ""), intern=FALSE, wait=TRUE)
				}else{
					suppressWarnings(success <- system(paste(shQuote(jags), if(!batch.jags) " <", " sim.", s, "/script.cmd > sim.", s, "/jagsoutput.txt 2>&1", sep=""), intern=TRUE, wait=TRUE))
				}
				
				retval <- "An error occured transferring the model output from the snow cluster"
				codafound <- FALSE
				# Wait for up to 5 secs for the coda files to be created:
				for(i in 1:5){
					Sys.sleep(1)
					if(file.exists("CODAindex.txt")) break
				}
				
				# Grab any new files from inside the simulation folder (there shouldn't be any folders):
				newfilenames <- list.files(paste("sim",s,sep="."))[! list.files(paste("sim",s,sep=".")) %in% names(simfiles)]
				newfiles <- vector('list',length=length(newfilenames))
				names(newfiles) <- newfilenames
				for(f in 1:length(newfiles)) newfiles[f] <- paste(readLines(paste("sim.",s,"/",newfilenames[f],sep=""), warn=FALSE), '\n', collapse="\n")
				
				if(cleanup){
					unlink(paste("sim",s,sep="."),recursive=TRUE)
					for(f in 1:length(mdfiles)){
						unlink(names(mdfiles)[f])
					}					
				}
				})
			
			return(list(retval=retval, newfiles=newfiles))
			
		}
		
		outfile <- ''
		if(identical(cl, NA)){
			outfile <- tempfile()
			cl <- parallel::makeCluster(n.sims, type=if(.Platform$OS.type=='unix') 'FORK' else 'PSOCK', outfile=outfile)
			on.exit({
				stopCluster(cl)
				unlink(outfile)
			})
		}
		
		clname <- paste(summary(cl)[1,'Class'], capture.output(print(cl)))
		clname <- gsub('forknode socket', 'Fork', clname)
		clname <- gsub('SOCKnode socket', 'PSOCK', clname)
		if(n.sims==1) swcat("Starting the simulation using a ", clname, "\n") else swcat("Starting ", n.sims, " simulations using a ", clname, "\n")		
		
		if(runjags.getOption('debug')==23){
			# Force recreation of the files to test it works OK:
			unlink(paste('sim',1:n.sims,sep='.'), recursive=TRUE)
			cat('Forcing file recreation for snow method\n')
		}
		
		success <- try({
			if(runjags.getOption('debug')>=100){
				warning('Using lapply in debug mode for snow')
				returns <- lapply(1:n.sims, f, files=files, jags=jags, batch.jags=batch.jags, silent.jags=silent.jags, path=thed)
			}else{
				returns <- parLapply(cl, 1:n.sims, f, files=files, jags=jags, batch.jags=batch.jags, silent.jags=silent.jags, path=thed)
			}
		}, silent=TRUE)
		if(inherits(success, 'try-error')){
			if(!file.exists(outfile))
				err <- paste("One or more snow sessions failed with the following error:\n",as.character(success),'\n',sep='')
			else
				err <- paste("One or more snow processes failed with the following error:\n",as.character(success),'\nThe worker log file (which may also help with debugging) is:\n',paste(readLines(outfile, warn=FALSE),collapse='\n'),sep='')
			
	    	stop(err,call.=FALSE)
		}
		
		charclass <- sapply(returns, function(x) return(length(x)==1 && class(x)=="character"))
		if(any(charclass)){
			cat("The following error(s) occured on one or more snow nodes:",unlist(returns[charclass]),sep="\n","")
		}
			
		retval <- "An error occured whilst writing new files to disk"
		
		success <- try({
			for(s in 1:n.sims){
				cwd <- getwd()
				if(!file.exists(paste("sim",s,sep="."))) dir.create(paste("sim",s,sep="."))
				setwd(paste("sim",s,sep="."))
				newfilenames <- names(returns[[s]]$newfiles)[! names(returns[[s]]$newfiles) %in% list.files()]
				# print(paste(length(newfilenames), "new files found for simulation", s))
				if(length(newfilenames)>0){
					for(f in 1:length(newfilenames)){
						cat(unlist(returns[[s]]$newfiles[newfilenames[f]]), file=newfilenames[f])
						# Close to make sure it's written out:
						close(file(newfilenames[f]))
					}
				}
				setwd(cwd)
			}
			})
		
		if(class(success)!="try-error"){
			retval <- TRUE
		}
		
	})
	
	return(retval)
	
}




runjags.background <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, rjags){
	
	if(n.sims==1) retval <- 'An unknown error occured while calling JAGS using the background method' else retval <- 'An unknown error occured while calling JAGS using the parallel background method'
	
	thed <- getwd()
	
  	bg.alert <- runjags.getOption('bg.alert')
  	
	if(runjags.getOption('bg.alert')!=""){
	  if (os == "windows"){
      lsfile <- 'lastsimlauncher.bat'
      cat(paste(shQuote(jags), if(!batch.jags) " <", " sim.", n.sims, "/script.cmd > sim.", n.sims, "/jagsoutput.txt 2>&1", sep=""),'\n',file=lsfile)      
      if(bg.alert=="beep"){
          bg.alert <- paste(shQuote(system.file('xgrid', 'beep.bat', package='runjags')),'\n',sep='')
		  #  Doesn't seem to produce 2 beeps....
#          bg.alert <- paste(shQuote(system.file('xgrid', 'beep.bat', package='runjags')),'\n',shQuote(system.file('xgrid', 'beep.bat', package='runjags')),sep='')
      }
      cat(bg.alert,'\n',file=lsfile,append=TRUE)
      
	  }else{
	    lsfile <- 'lastsimlauncher.sh'
	    cat(paste(shQuote(jags), if(!batch.jags) " <", " sim.", n.sims, "/script.cmd > sim.", n.sims, "/jagsoutput.txt 2>&1", sep=""),'\n',file=lsfile)      
	    if(bg.alert=="beep"){
	      bg.alert <- "printf '\\007'\nsleep 1\nprintf '\\007'\n"
	    }
	    cat(bg.alert,'\nexit 0\n',file=lsfile,append=TRUE)
	    Sys.chmod(lsfile)
	  }
	}	  
	  
	tryCatch({

		if(n.sims==1) swcat("Starting the simulation in the background...\n") else swcat("Starting ", n.sims, " simulations in the background...\n")
		
		success <- numeric(n.sims)
		for(s in 1:n.sims){
		  	command <- paste(shQuote(jags), if(!batch.jags) " <", " sim.", s, "/script.cmd > sim.", s, "/jagsoutput.txt 2>&1", sep="")      
      		if(s == n.sims){
        		command <- shQuote(file.path(getwd(),lsfile))
      	  	}
			
			if(runjags.getOption('debug'))
				swcat('Simulation ', s, ' of ', n.sims, ' - ', os, ' command:  ', command, '\n', sep='')	
			
			if(os == "windows"){
				success[s] <- shell(command, intern=FALSE, wait=FALSE)
			}else{
				success[s] <- system(command, intern=FALSE, wait=FALSE)
			}
		}
		
		Sys.sleep(1)
	
		retval <- FALSE
		
		if(n.sims==1) swcat("The JAGS process is now running in the background\n") else swcat("The JAGS processes are now running in the background\n")
	})
	
	return(retval)
}


runjags.interruptible <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, rjags){
	
	if(runjags.getOption('debug'))
		swcat('Using the following simulation directory: ', getwd(), '\n')
	
	swcat("Calling the simulation...\n")

	retval <- 'An unknown error occured while calling JAGS using the interruptible method'
	
	stopifnot(n.sims==1)
	
	#filecon <- file('scriptlauncher.sh', 'w')
	cat('#!/bin/sh
	i=$1
	', shQuote(jags), if(!batch.jags) ' <', ' sim.$i/script.cmd > sim.$i/jagsoutput.txt 2>&1 &

	echo $! > sim.$i/jagspid.txt
	', sep='', file='scriptlauncher.sh')
	#close(filecon)

	Sys.chmod('scriptlauncher.sh')

	thed <- getwd()
	
	tryCatch({

		
		if(os == "windows"){
		
			tasks <- system('TASKLIST', intern=TRUE)

			allpid=as.numeric(unlist(lapply(tasks[grepl("jags-terminal", tasks) | grepl("JAGS-T~1.EXE", tasks) ], function(x){
				chars <- strsplit(x, split=" ")[[1]]
				return(chars[chars!=""][2])
				})))

			output <- shell(paste(shQuote(jags), if(!batch.jags) " <", " sim.1/script.cmd > sim.1/jagsoutput.txt 2>&1", sep = ""), intern=FALSE, wait=FALSE)

			tasks <- system('TASKLIST', intern=TRUE)
			newpid=as.numeric(unlist(lapply(tasks[grepl("jags-terminal", tasks) | grepl("JAGS-T~1.EXE", tasks) ], function(x){
				chars <- strsplit(x, split=" ")[[1]]
				return(chars[chars!=""][2])
				})))

			pid <- newpid[! newpid %in% allpid]
			
		}else{
			
			success <- system('./scriptlauncher.sh 1', wait=TRUE, intern=FALSE)	

			# Allow simulation to start before we look for the PID:
			Sys.sleep(0.1)
			if(!file.exists('sim.1/jagspid.txt')) Sys.sleep(1)
			if(!file.exists('sim.1/jagspid.txt')) Sys.sleep(1)
			
			suppressWarnings(output <- readLines('sim.1/jagspid.txt', warn=FALSE))
			pid <- output[1]
			if(as.numeric(pid)!=as.integer(pid) | pid=="") stop("A problem occured when reading the output of a started process")
	
		}
		
		output <- tailf('sim.1/jagsoutput.txt', refresh=jags.refresh, start=1, min.static=2, stop.text=getstoptexts(), print=!silent.jags, return=TRUE)
	
		if(output$interrupt){
			
			retval  <- "The JAGS process was terminated by the user"
			
			if(os == 'windows'){
				# Make sure length(pid) == 1
				if(length(pid)!=1){
					warning("Unable to identify correct JAGS process to terminate; no processes have been killed (the JAGS model will continue until it has finished or you restart your computer)", call.=FALSE)
				}else{
					suppressWarnings(killout <- system(paste('taskkill /F /PID ', pid, sep=''), intern=TRUE))
				}

			}else{
				system(paste('kill ', pid, sep=''), ignore.stderr=TRUE)
			}
			
		}else{
			retval <- TRUE	
		}
	})
	
	flush.console()
	swcat("\n")
	setwd(thed)
	return(retval)	

}

runjags.parallel <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, rjags){
	
	if(n.sims==1) swcat("Calling the simulation using the parallel method...\n") else swcat("Calling ", n.sims, " simulations using the parallel method...\n")

	retval <- 'An unknown error occured while calling JAGS using the parallel method'
	
	#filecon <- file('scriptlauncher.sh', 'w')
	cat('#!/bin/sh
	i=$1
	', shQuote(jags), if(!batch.jags) ' <', ' sim.$i/script.cmd > sim.$i/jagsoutput.txt 2>&1 &

	echo $! > sim.$i/jagspid.txt
	', sep='', file='scriptlauncher.sh')
	#close(filecon)

	Sys.chmod('scriptlauncher.sh')

	thed <- getwd()
	
	
	tryCatch({

		if(os == "windows"){
		
			tasks <- system('TASKLIST', intern=TRUE)

			allpid=as.numeric(unlist(lapply(tasks[grepl("jags-terminal", tasks) | grepl("JAGS-T~1.EXE", tasks) ], function(x){
				chars <- strsplit(x, split=" ")[[1]]
				return(chars[chars!=""][2])
			})))

			for(s in 1:n.sims){
 				output <- shell(paste(shQuote(jags), if(!batch.jags) " <", " sim.", s, "/script.cmd > sim.", s, "/jagsoutput.txt 2>&1", sep = ""), intern=FALSE, wait=FALSE)
				# Allow simulation to start before we move onto the next one:
				Sys.sleep(0.1)
			}
						
			# Allow first simulation to start before we look for the PID:
			tries <- 0
			repeat{
				if(file.exists('sim.1/jagsoutput.txt')) break
				Sys.sleep(0.5)
				tries <- tries +1
				if(tries==6) stop("Timed out waiting for output to appear from JAGS")
			}
			
			tasks <- system('TASKLIST', intern=TRUE)

			newpid=as.numeric(unlist(lapply(tasks[grepl("jags-terminal", tasks) | grepl("JAGS-T~1.EXE", tasks) ], function(x){
				chars <- strsplit(x, split=" ")[[1]]
				return(chars[chars!=""][2])
			})))

			pid <- newpid[! newpid %in% allpid]
		
		}else{
			pid <- character(n.sims)
			# Start with last simulation first, so the first simulation will be the last to finish most of the time:
			for(s in n.sims:1){
				success <- suppressWarnings(system(paste('./scriptlauncher.sh ', s, sep=''), wait=TRUE, intern=FALSE))
				# Allow simulation to start before we look for the PID:
				tries <- 0
				repeat{
					if(file.exists(paste('sim.', s, '/jagspid.txt', sep=''))) break
					Sys.sleep(0.5)
					tries <- tries +1
					if(tries==6) stop("Timed out waiting for output to appear from JAGS")
				}

				suppressWarnings(output <- readLines(paste('sim.', s, '/jagspid.txt', sep=''), warn=FALSE))
				pid[s] <- output[1]
				if(as.numeric(pid[s])!=as.integer(pid[s]) | pid[s]=="") stop("A problem occured when reading the output of a started process")		
			}
		}

		# Loop through and follow chains sequentially until they've all finished:
		s <- 1
		repeat{
			if(!silent.jags) swcat("Following the progress of chain ", s, " (the program will wait for all chains to finish before continuing):\n", sep="")			
			output <- tailf(paste('sim.', s, '/jagsoutput.txt', sep=''), refresh=jags.refresh, start=1, min.static=2, stop.text=getstoptexts(), print=!silent.jags, return=TRUE)
			if(!silent.jags) swcat("\n")
			if(output$interrupt) break

			# The first sim should be finished last ... but give the others up to 5 seconds to catch up before switching to them...
			tries <- 0
			repeat{
				# Check to see which chains have (a) started, and (b) finished:
				simsdone <- sapply(1:n.sims, function(x) return(file.exists(paste('sim.', x, '/jagsoutput.txt', sep='')) && 
					any(sapply(getstoptexts(), grepl, x=paste(readLines(paste('sim.', x, '/jagsoutput.txt', sep=''),warn=FALSE),collapse="\n")))))
				
				if(all(simsdone)) break
				tries <- tries +1
				if(tries==5) break
				Sys.sleep(1)
			}

			if(all(simsdone)){
				swcat("All chains have finished\n")
				break
			}

			# If one still hasn't finished, start watching it:
			s <- which(!simsdone)[1]
		}					
				
		if(output$interrupt){

			retval  <- "The JAGS process was terminated by the user"
			if(os == 'windows'){
				# Make sure length(pid) == n.sims
				if(length(pid)!=n.sims){
					warning("Unable to identify correct JAGS processes to terminate; no processes have been killed (the JAGS model will continue until it has finished or you restart your computer)", call.=FALSE)
				}else{
					for(k in pid){
						suppressWarnings(killout <- system(paste('taskkill /F /PID ', k, sep=''), intern=TRUE))
					}
				}
				
			}else{
				for(s in 1:n.sims){
					system(paste('kill ', pid[s], sep=''), ignore.stderr=TRUE)
				}
			}			

		}else{
			
			retval <- TRUE
				
		}
	
	})
	
	flush.console()	
	return(retval)

}


runjags.rjags <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, rjags){
	
	# rjags object should have been compiled by extend.jags or autoextend.jags if this was necessary - and the modules/factories loaded just before that
	extra.options <- rjags[names(rjags)!='rjags']
	rjags <- rjags$rjags

	if(silent.jags) extra.options$progress.bar <- "none"
	
	if(!silent.jags) swcat("Calling the simulation using the rjags method...\n")

	# We may need this current state again later:
	if(runjags.getOption('repeatable.methods')) inits <- rjags$state(internal=TRUE)
	
	needsadapting <- checkadaptrequired(rjags)
	# Model is guaranteed to be compiled, but as.jags didn't adapt it - so need to check it is adapted (and/or do adapting):
	if(extra.options$adapt>0){
		# After the model is compiled, adapt(0) will either return an error or TRUE
		if(needsadapting){	
			if(!silent.jags)
				swcat("  Adapting the model for ", format(extra.options$adapt,scientific=FALSE), " iterations...\n",sep="")
		
			flush.console()
			by <- if(is.na(extra.options$by)) min(100, extra.options$adapt/50) else extra.options$by
			adaptdone <- rjags::adapt(rjags,n.iter=extra.options$adapt,end.adaptation=TRUE,progress.bar=extra.options$progress.bar,by=by)
			flush.console()
		
			if(!adaptdone){
				if(runjags.getOption('adapt.incomplete')=='error')
					stop(paste('The adaptation phase of the model was not completed in ', extra.options$adapt, ' iterations - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
				if(runjags.getOption('adapt.incomplete')=='warning')
					warning(paste('The adaptation phase of the model was not completed in ', extra.options$adapt, ' iterations, so the current samples may not be optimal - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
			}
		}else{
			if(!silent.jags && runjags.getOption('adapt.incomplete')!='silent')
				swcat('Note: the model did not require adaptation\n')
		}		
	}else{
		if(needsadapting && extra.options$burnin <1000){   # If burnin is long enough, assume that adaptation is OK until I can do adapt(,0)
			if(runjags.getOption('adapt.incomplete')=='error')
				stop(paste('The model requires an adaptation phase and the burnin period is not long enough to guarantee adaptation - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
			if(runjags.getOption('adapt.incomplete')=='warning')
				warning(paste('The model requires an adaptation phase and the burnin period is not long enough to guarantee adaptation, so the current samples may not be optimal - try increasing the number of iterations to the "adapt" argument', sep=''), call.=FALSE)
		}
	}	

	if(extra.options$burnin>0){
		if(!silent.jags)
			swcat("  Burning in the model for ", format(extra.options$burnin,scientific=FALSE), " iterations...\n",sep="")
		flush.console()
		by <- if(is.na(extra.options$by)) min(100, extra.options$burnin/50) else extra.options$by
		update(rjags,n.iter=extra.options$burnin,progress.bar=extra.options$progress.bar,by=by)
		flush.console()
	}
	
	# The behaviour of rjags method for dic stats is slightly different to the others, but should appear to be the same:
	# If just DIC is specified, do like I have always done BUT don't return pd and deviance unless specifically requested
	# If ped is specified, call dic.samples
	# If full.pd is specified, call dic.samples
	
	# Remove full.pd and popt monitors
	monitor <- extra.options$monitor[! extra.options$monitor %in% c("popt","full.pd","dic","ped")]
	monitor[monitor=="pd"] <- "pD"
	if(any(extra.options$monitor=='dic'))
		monitor <- unique(c(monitor, 'pD', 'deviance'))	
	
	by <- if(is.na(extra.options$by)) min(100, extra.options$burnin/50) else extra.options$by
	if(!silent.jags) swcat("  Running the model for ", format(extra.options$sample,scientific=FALSE), " iterations", if(runjags.getOption('debug')>5) paste(" (updating by ", by, ")", sep=""), "...\n",sep="")
	flush.console()	
	samples <- rjags::jags.samples(rjags,variable.names=monitor,n.iter=extra.options$sample,progress.bar=extra.options$progress.bar, thin=extra.options$thin,by=by)
	
	# Do this BEFORE our dummy call of 1 iteration
	end.state <- sapply(rjags$state(internal=TRUE),dump.format)	
	
	# This is just a dummy call so that we can get the names of the variables:
	suppressWarnings(varnames <- varnames(rjags::coda.samples(rjags,variable.names=monitor[monitor!="pD"],n.iter=2,progress.bar="none", thin=1)))
	flush.console()
	mcmcout <- lapply(samples[names(samples)!='pD'], as.mcmc.list)
	
	nvar <- length(varnames)
  
  	niter <- sapply(mcmcout,niter)
	if(!all(niter==niter[1])) stop("An error occured with the rjags method - variables returned with differing numbers of iterations")
	niter <- niter[1]

	########
	# For deviance.sum:
	sum.pd <- NA
	sum.popt <- NA
	sum.deviance <- NA
	# For deviance.table:
	pd <- NA
	popt <- NA
	deviance <- NA
	# For $pd:
	fullsumpd <- NA
	########
	
	# rjags (verison 3) returns variables in alphabetical order regardless of input order, so re-order here:
	reindex <- matchvars(monitor[monitor!="pD"],varnames,exactneeded=TRUE,testfound=FALSE)
	varnames <- varnames[reindex]
	
	mcmc <- vector('list',length=extra.options$n.chains)
	for(i in 1:extra.options$n.chains){
		mcmc[[i]] <- mcmc(do.call('cbind',lapply(mcmcout, function(x) return(x[[i]])))[,reindex,drop=FALSE], start=extra.options$burnin+1, thin=extra.options$thin)
		dimnames(mcmc[[i]]) <- list(1:niter, varnames)
	}
	if('deviance'%in%varnames){
	  sum.deviance <- mean(sapply(1:extra.options$n.chains, function(x) return(mean(mcmc[[x]][,'deviance']))))
	}
	if(!'deviance'%in%extra.options$monitor && any(varnames=='deviance')){
		for(i in 1:extra.options$n.chains){
		  mcmc[[i]] <- mcmc[[i]][,-which(varnames=='deviance'),drop=FALSE]
		}
	}
	mcmc <- as.mcmc.list(mcmc)
		
	fullsumpd <- NA
	if(any(c("pD",'full.pd')%in%names(samples))){
		fullsumpd <- mcmc(matrix(unlist(samples[names(samples)=='pD']),ncol=1,dimnames=list(NULL,"pD")), start=extra.options$burnin+1, thin=extra.options$thin)
		dimnames(fullsumpd) <- list(1:niter(fullsumpd), dimnames(fullsumpd)[[2]])
		
		# We want the mean of the sum pD monitor here:
		sum.pd <- mean(fullsumpd)
		if(!'full.pd'%in%extra.options$monitor)
			fullsumpd <- NA
	}
	
	if(any(c('ped','popt')%in%extra.options$monitor)){
		if(runjags.getOption('repeatable.methods')){
		  	if(!silent.jags) swcat("  Re-compiling the model with original starting values...\n",sep="")
			
			# Retrieve model, data and initial values with RNG states:
			model <- paste(rjags$model(),collapse="\n")
			data <- rjags$data()
			# Prevents problems with data/inits reloading:
			data <- data[names(data) %in% extra.options$origdatanames]
			
			n.chains <- length(inits)
			tmodel <- textConnection(model)

			rjags <- rjags::jags.model(file=tmodel, data=data, inits=inits, n.chains, n.adapt=0, quiet=TRUE)
			close(tmodel)
		
			if(extra.options$adapt>0 && needsadapting){
				if(!silent.jags) swcat("  Re-adapting the model for ", format(extra.options$adapt,scientific=FALSE), " iterations...\n",sep="")
				flush.console()
				by <- if(is.na(extra.options$by)) min(100, extra.options$adapt/50) else extra.options$by				
				adaptdone <- rjags::adapt(rjags,n.iter=extra.options$adapt,end.adaptation=TRUE,progress.bar=extra.options$progress.bar, by=by)
			}
		
			if(extra.options$burnin>0){
				if(!silent.jags) swcat("  Re-burning in the model for ", format(extra.options$burnin,scientific=FALSE), " iterations...\n",sep="")
				flush.console()
				by <- if(is.na(extra.options$by)) min(100, extra.options$burnin/50) else extra.options$by
				update(rjags,n.iter=extra.options$burnin,progress.bar=extra.options$progress.bar, by=by)
				flush.console()
			}

			flush.console()
		}
		
	  	if(!silent.jags) swcat("  Extending ", format(extra.options$sample,scientific=FALSE), " iterations for pOpt/PED estimates...\n",sep="")
	  	flush.console()
		by <- if(is.na(extra.options$by)) min(100, extra.options$sample/50) else extra.options$by
		
	  	dics <- rjags::dic.samples(rjags, n.iter=extra.options$sample, thin=extra.options$thin, type='popt', progress.bar=extra.options$progress.bar, by=by)
	  	flush.console()	  
    	deviance <- as.matrix(dics$deviance)
    	popt <- as.matrix(dics$penalty)	  
		dn <- dimnames(as.matrix(dics$deviance))[[1]]

		# Ensures mean pd is not replicated:
		sum.deviance <- sum(deviance)
		sum.popt <- sum(popt)
	}
	if(any(c('pd')%in%extra.options$monitor) || (runjags.getOption('repeatable.methods') && 'dic'%in%extra.options$monitor)){
		if(runjags.getOption('repeatable.methods')){
		  	if(!silent.jags) swcat("  Re-compiling the model with original starting values...\n",sep="")
			
			# Retrieve model, data and initial values with RNG states:
			model <- paste(rjags$model(),collapse="\n")
			data <- rjags$data()
			# Prevents problems with data/inits reloading:
			data <- data[names(data) %in% extra.options$origdatanames]
			
			n.chains <- length(inits)
			tmodel <- textConnection(model)

			rjags <- rjags::jags.model(file=tmodel, data=data, inits=inits, n.chains, n.adapt=0, quiet=TRUE)
			close(tmodel)
		
			if(extra.options$adapt>0 && needsadapting){
				if(!silent.jags) swcat("  Re-adapting the model for ", format(extra.options$adapt,scientific=FALSE), " iterations...\n",sep="")
				flush.console()
				by <- if(is.na(extra.options$by)) min(100, extra.options$adapt/50) else extra.options$by				
				adaptdone <- rjags::adapt(rjags,n.iter=extra.options$adapt,end.adaptation=TRUE,progress.bar=extra.options$progress.bar,by=by)
			}
		
			if(extra.options$burnin>0){
				if(!silent.jags) swcat("  Re-burning in the model for ", format(extra.options$burnin,scientific=FALSE), " iterations...\n",sep="")
				flush.console()
				by <- if(is.na(extra.options$by)) min(100, extra.options$burnin/50) else extra.options$by
				update(rjags,n.iter=extra.options$burnin,progress.bar=extra.options$progress.bar, by=by)
				flush.console()
			}

			flush.console()
		}

	  	if(!silent.jags) swcat("  Extending ", format(extra.options$sample,scientific=FALSE), " iterations for pD/DIC estimates...\n",sep="")
	  	flush.console()
		by <- if(is.na(extra.options$by)) min(100, extra.options$sample/50) else extra.options$by
	  	dics <- rjags::dic.samples(rjags, n.iter=extra.options$sample, thin=extra.options$thin, type='pD', progress.bar=extra.options$progress.bar, by=by)
	  	flush.console()	  
    	deviance <- as.matrix(dics$deviance)
    	pd <- as.matrix(dics$penalty)	  
		dn <- dimnames(as.matrix(dics$deviance))[[1]]

		sum.deviance <- sum(deviance)
		sum.pd <- sum(pd)
	}
	
	deviance.sum=deviance.table <- NA
	if(any(!is.na(c(deviance, pd, popt)))){
		deviance.table <- cbind(deviance, pd, popt)
		dimnames(deviance.table) <- list(dn,c('mean.deviance','mean.pD','mean.pOpt'))
	}
	if(any(!is.na(c(sum.deviance, sum.pd, sum.popt)))){
		deviance.sum <- c(sum.deviance, sum.pd, sum.popt)
		names(deviance.sum) <- c('sum.mean.deviance', 'sum.mean.pD', 'sum.mean.pOpt')
	}
	
	success <- try({
		samplers <- rjags::list.samplers(rjags)
		samplers <- as.data.frame(t(matrix(unlist(lapply(1:length(samplers), function(x) return(rbind(Index.No=x, Sampler=names(samplers)[x], Node=samplers[[x]])))),nrow=3)), stringsAsFactors=FALSE)
		dimnames(samplers)[[2]] <- c('Index.No', 'Sampler', 'Node')
		samplers$Index.No <- as.numeric(samplers$Index.No)	
		}, silent=TRUE)
	if(inherits(success, 'try-error')){
		# warning('There was an error retrieving the samplers',call.=FALSE)
		# This is dealt with elsewhere
		samplers <- NA
	}
	
	if(!silent.jags) swcat("Simulation complete\n")
	
	return(list(complete=TRUE, mcmc=mcmc, deviance.table=deviance.table, deviance.sum=deviance.sum, pd=fullsumpd, end.state=end.state, samplers=samplers))

}

runjags.xgrid <- function(jags, silent.jags, jags.refresh, batch.jags, os, libpaths, n.sims, jobname, cl, remote.jags, command, customart, jagspath, submitandstop, max.threads, mgridpath, hostname, password){
	
	swcat("Submitting the simulation to Xgrid... (this may take some time)\n")
		
	retval <- 'An unknown error occured while calling JAGS using the xgrid method'
	
	# max.threads is only passed through for consistency, and should be redundant by n.sims:
	stopifnot(n.sims<=max.threads)
			
	#filecon <- file('scriptlauncher.sh', 'w')
	cat('#!/bin/sh
	i=$1
	', shQuote(jags), if(!batch.jags) ' <', ' sim.$i/script.cmd > sim.$i/jagsoutput.txt 2>&1 &

	echo $! > sim.$i/jagspid.txt
	', sep='', file='scriptlauncher.sh')
	#close(filecon)

	Sys.chmod('scriptlauncher.sh')
	
	thed <- getwd()
	
	tryCatch({


		cat(customart, file="customart.sh")

	# cd doesn't work on some nodes, so the working directory is left where it is now and paths adjusted accordingly in the script
	#	', if(method.options$separate.tasks) 'cd sim.$1', '
		cat('#!/bin/sh


		pid=$$

		', if(!silent.jags) 'echo "" > sim.$1/jagsout.txt', ' 
		', if(!silent.jags) 'echo "" > sim.$1/jagsout.txt', ' 
		', if(n.sims>1) '( ( echo "Chain "$1":" 2>&1 1>&3 | tee -a sim.$1/jagserror.txt) 3>&1 1>&2) 2>&1 | tee -a sim.$1/jagsout.txt', if(n.sims>1 && silent.jags) '', '', '

		( ( (', jagspath, ' < sim.$1/script.cmd; echo $! > sim.$1/.retstat.$pid) 2>&1 1>&3 | tee -a sim.$1/jagserror.txt) 3>&1 1>&2) 2>&1 | tee -a sim.$1/jagsout.txt', if(silent.jags) '', '

		', if(n.sims>1) '( ( echo "" 2>&1 1>&3 | tee -a sim.$1/jagserror.txt) 3>&1 1>&2) 2>&1 | tee -a sim.$1/jagsout.txt', if(n.sims>1 & silent.jags) '', '', '

		# This makes sure the process has finished before continuing:
		wait
	
		returnstat=`cat < sim.$1/.retstat.$$`
		rm sim.$1/.retstat.$$
	
		', if(!silent.jags) 'echo ""', '
	
		exit $returnstat
			', sep='', file=jobname)
			#close(filecon)
	
		Sys.chmod(jobname)

		
		cat('#!/bin/sh		
	cmd="', jobname, '" 
	ntasks=', n.sims, '
	job=1
	indir="', thed, '"
	# Start the process in the background so we can get the pid:
	', command, ' &
	# Echo the pid to a file - we might need to kill it if interrupted:
	echo $! > mgridpid.txt
	# And wait for mgrid to finish:
	wait
	rm mgridpid.txt
	exit 0
		', sep='', file='starter.sh')
	#}

		Sys.chmod('starter.sh')
		
		# For some reason wrapping the actual system call inside a shell script fixes the display issues in Aqua, so no need to use tailf workaroud here:
		success <- system('( ./starter.sh 2>&3 | tee .starterout.txt) 3>&2 | tee starteroutput.txt', intern=FALSE, wait=TRUE)

		# This should be remvoed by the shell script, if it's still there the shell script didn't finish:
		if(file.exists("mgridpid.txt")){
			interrupt <- TRUE
			pid <- as.numeric(readLines('mgridpid.txt', warn=FALSE))
			success <- system(paste("kill ", pid, sep=""))
			return("The process was terminated while submitting or waiting for the results of the xgrid job")
		}else{
				
			if(success!=0){
				return("An unknown error occured while starting the xgrid command")
			}
	
			if(file.exists('jobid.txt')) joboutput <- paste(readLines('jobid.txt', warn=FALSE), collapse='\n') else joboutput <- paste(readLines('starteroutput.txt', warn=FALSE), collapse='\n')
				
			results <- list()
			
			if(submitandstop){
				jobnum <- gsub('[^[:digit:]]', '', paste(joboutput, collapse=''))
				if(jobnum=='' | as.numeric(jobnum)>10^6){
			#				cat(readLines("starteroutput.txt", warn=FALSE))
					return("There was an error submitting your job to Xgrid - no Job ID was returned")			
				}
				swcat('Your job (name: "', jobname, '";  ID: "', jobnum, '") has been succesfully uploaded to xgrid\n', sep='')
				results$complete <- FALSE
				results$jobid <- jobnum
			}else{
				results$complete <- TRUE
				results$jobid <- NA
			}
			
			retval <- results
	
		}
	
	})
	
	return(retval)
	
}

runjags.start <- function(model, monitor, data, inits, modules, factories, burnin, sample, adapt, thin, tempdir, dirname, method, method.options, internal.options=list()){
	
	# Reset failedjags stuff:
	failedjags$model <- NA
	failedjags$data <- NA
	failedjags$inits <- NA
	failedjags$output <- NA
	failedjags$end.state <- NA
	
	modules <- checkmodfact(modules, 'module')
	factories <- checkmodfact(factories, 'factory')
	
	# In case noread.monitor and monitor partly overlap:
	monitor <- unique(monitor)
	doublecounted <- function(x) return(grepl('\\[.*\\]',x) && any(gsub('\\[.*\\]', '', x)==monitor))
	monitor <- monitor[!sapply(monitor,doublecounted)]
	 
	chains=n.chains <- length(inits)

	updates <- sample
	if(sample<1) stop("The specified value for sample must be a positive integer")

	# Method for calling JAGS is now encapsulated within method=function() and method.options=list()
	# jags, silent.jags, tempdir, method, jags.refresh, batch.jags){

	# Some methods don't require writing files:
	writefiles <- TRUE
	
	# Checking of DIC and JAGS availability etc will already have been done, but we need to check vailability of methods (and xgrid options)
	if(class(method)!="function"){
		
		method <- getrunjagsmethod(method)
		
		##### Internal.options is ONLY used by Xgrid - get rid of it once Xgrid is gone
		if(method=="xgrid"){			
			xgrid.options <- method.options
			method.options <- internal.options
			cl = inputsims = remote.jags = rjags = by = progress.bar <- NA			
		}else{
			requirednames <- c("n.sims","cl","remote.jags","by","progress.bar","jags","silent.jags","jags.refresh","batch.jags")
			# The method.options are now being checked and handled by extend.jags - this should always be true (note rjags is not required as we remove it occasionally)
			if(!all(requirednames%in%names(method.options)) || any(is.null(method.options[which(names(method.options)%in%requirednames)]))){
				if(runjags.getOption('debug')){
					cat('**********  Method options dont match\n')
					browser()
				}
				warning('There was an unexpected problem with the method options list', call.=FALSE)
				method.options <- method.options[names(method.options)%in%c("n.sims","cl","remote.jags","rjags","by","progress.bar","jags","silent.jags","jags.refresh","batch.jags")]
			}

			# This is a bit clunky at the moment - what was method.options (the 5 below) are handled differently to what was input.options (now in method.options, but previously method.options was reset to input.options):
			inputsims <- method.options$n.sims
			cl <- method.options$cl 
			remote.jags <- method.options$remote.jags 
			rjags <- method.options$rjags 
			by <- method.options$by 
			progress.bar <- method.options$progress.bar
						
			# Need to remoe these as they are re-appended later:
			method.options <- method.options[!names(method.options)%in%c('n.sims','cl','remote.jags','rjags','by','progress.bar')]
		}
		
		# Use the cluster to get available cores otherwise parallel:
		if(any(grepl('cluster',class(cl)))){
			ncores <- length(cl)
		}else{
			ncores <- max(2, suppressWarnings(parallel::detectCores()))  # never less than 2 cores!
		}
		
		#####  The OLD code for comparison:
#		cl = inputsims = remote.jags = rjags = by = progress.bar <- NA
#		if(method=="xgrid"){			
#			xgrid.options <- method.options
#			method.options <- internal.options
#		}else{
#			if(!is.null(method.options$n.sims) && is.numeric(method.options$n.sims)) inputsims <- method.options$n.sims
#			if(!is.null(method.options$cl) && any(grepl('cluster',class(method.options$cl)))) cl <- method.options$cl 
#			if(!is.null(method.options$remote.jags)) remote.jags <- method.options$remote.jags 
#			if(!is.null(method.options$rjags)) rjags <- method.options$rjags 
#			if(!is.null(method.options$by)) by <- method.options$by 
#			if(!is.null(method.options$progress.bar)){
#				progress.bar <- method.options$progress.bar
#			}else{
#				progress.bar <- unlist(options('jags.pb'))
#			}								
#			if(!identical(method.options, list()) && !identical(names(method.options), c('jobname', 'command', 'customart', 'jagspath', 'submitandstop', 'max.threads', 'mgridpath', 'hostname', 'password')) && !all(names(method.options)%in%c("n.sims","cl","remote.jags","rjags","by","progress.bar"))) warning("The supplied value for method.options is ignored when using an inbuilt method (except for 'n.sims' which can be provided for parallel methods, 'cl' and 'remote.jags' which can be provided for the snow method, and 'by' and 'progress.bar' which can be supplied for the rjags method)")
#			method.options <- internal.options
#		}
		
		
		
		
		# Ignore tempdir argument if using rjagsmethod
		if(method%in%runjagsprivate$rjagsmethod) tempdir <- TRUE
		
		# Getting a strange no visible binding error from R CMD check - can't track it down but this should cure it:
		jags <- method.options$jags
		
		if((method=='parallel' | method=='interruptible') &.Platform$OS.type!='windows' ) if(length(suppressWarnings(system('/bin/ps', intern=TRUE, ignore.stderr=TRUE)))==0 & Sys.info()['user']=='nobody'){
			warning('The interruptible and parallel methods are unavailable when running over xgrid, using the simple method instead')
			method <- 'simple'
		}
		if((method=='parallel' | method=='interruptible') &.Platform$OS.type=='windows'){
			if(Sys.which('TASKLIST')[1]==""){
				warning("Parallel and interruptible methods aren't available on your machine because the TASKLIST system command couldn't be found; switching to the simple method", call.=FALSE)
				method <- 'simple'
			}		
			if(Sys.which('TASKKILL')[1]==""){
				warning("Parallel and interruptible methods aren't available on your machine because the TASKKILL system command couldn't be found; switching to the simple method", call.=FALSE)
				method <- 'simple'
			}
			ret <- system2('TASKLIST',stdout=FALSE,stderr=FALSE)
			if(ret!=0){
				warning("Parallel and interruptible methods aren't available on your machine because testing the TASKLIST system command produced an error; switching to the simple method", call.=FALSE)
				method <- 'simple'			
			}
			ret <- system2('TASKKILL','/?',stdout=FALSE,stderr=FALSE)
			if(ret!=0){
				warning("Parallel and interruptible methods aren't available on your machine because testing the TASKKILL system command produced an error; switching to the simple method", call.=FALSE)
				method <- 'simple'			
			}
		}		
	
		jags.status <- testjags(method.options$jags, silent=TRUE)
		if(jags.status$JAGS.available==FALSE){
			if(jags.status$os=="windows"){
				# Try it again - sometimes this seems to clear it up:
				Sys.sleep(0.2)
				jags.status <- testjags(method.options$jags, silent=TRUE)
			}		
			jags <- jags.status$JAGS.path
		
			if(jags.status$JAGS.available==FALSE){			
				swcat("Unable to call JAGS using '", jags, "' - try specifying the path to the JAGS binary as the jags argument, or installing the rjags package.  Use the testjags() function for more detailed diagnostics.\n", sep="")
				stop("Unable to call JAGS", call.=FALSE)
			}
		}
		jags <- jags.status$JAGS.path
		
		if(!jags.status$JAGS.found && ! method%in% c("snow",runjagsprivate$rjagsmethod)){
			swcat("Unable to call JAGS using '", jags, "' - try specifying the path to the JAGS binary as the jags argument, or using the rjags method.  Use the testjags() function for more detailed diagnostics.\n", sep="")
			stop("Unable to call JAGS", call.=FALSE)
		}
		if(method%in%runjagsprivate$rjagsmethod){
			loadandcheckrjags(stop=TRUE)
		}
	
		if(jags.status$JAGS.major=="unknown" | is.na(jags.status$JAGS.major)){
			warning('Unable to verify the version number of JAGS.  If any functions do not work as expected, you could try checking your JAGS installation for problems.')
			jags.status$JAGS.major <- Inf
		}
	
		if(length(grep('base::Mersenne-Twister', inits)>0) & as.numeric(jags.status$JAGS.major) < 2) warning('Using the RNG "base::Mersenne-Twister" (used by default for chain 4) may cause problems with restarting subsequent simulations using the end state of previous simulations due to a bug in JAGS version 1.x.  If you encounter the error "Invalid .RNG.state", please update JAGS to version 2.x and try again.  Or, you can change the random number generator by changing the .RNG.name to (for example) "base::Super-Duper" and remove the .RNG.state element of the list.')
		
		if(tempdir){
			if(is.na(dirname)) dirname <- "runjagsfiles"
			temp.directory <- tempfile(dirname)
			jobname <-  if(method=='xgrid') paste("xgrid.",dirname,sep="") else temp.directory
			dir.create(temp.directory)
		}else{
			# Change jobname to match folder name:
			if(method=='xgrid'){
				jobname <- new_unique(xgrid.options$jobname, touch=TRUE, type='folder')
			}else{
				if(is.na(dirname)){
					dirname <- "runjagsfiles"
					jobname <- new_unique(dirname, touch=TRUE, type='folder')
				}else{
					jobname <- new_unique(dirname, touch=TRUE, type='folder')
					if(jobname!=dirname) warning(paste("The specified folder name was already taken - saving to '", jobname, "' instead", sep=""))
				}
			}
			if(jobname=="Directory not writable"){
				stop("Directory not writable", call.=FALSE)
			}
			temp.directory <- file.path(getwd(), jobname)
		}

		strmethod <- method
		extramethodargs <- list(sim.chains=list())
		if(strmethod=="simple"){
			n.sims <- 1
			method <- runjags.simple
		}
		if(strmethod=="interruptible"){
			n.sims <- 1
			method <- runjags.interruptible
		}
		if(strmethod=="parallel"){			
			if(is.na(ncores)){
				ncores <- 2
				warning("Unable to detect the available number of cores on your machine - using a maximum of 2 cores as a default")
			}
			n.sims <- min(n.chains, ncores)
			if(!identical(NA, inputsims)){
				if(inputsims!=as.integer(inputsims) || inputsims <1) warning("The supplied n.sims option is not a positive integer and was ignored") else n.sims <- min(n.chains, inputsims)
			}
			method <- runjags.parallel
		}
		if(strmethod=="snow"){
			if(is.na(ncores)){
				ncores <- 2
				warning("Unable to detect the available number of cores on your machine - using a maximum of 2 cores as a default")
			}
			n.sims <- min(n.chains, ncores)
			if(!is.na(inputsims)){
				if(inputsims!=as.integer(inputsims) || inputsims <1) warning("The supplied n.sims option is not a positive integer and was ignored") else n.sims <- min(n.chains, inputsims)
			}
			method <- runjags.snow
			rjt <- remote.jags
			if(class(rjt)=="function"){
				s <- try(rjt <- rjt())
				if(class(s)=="try-error") stop("The function supplied to remote.jags must take no arguments")
			}
			if(!is.character(rjt) && length(rjt)!=1)
				stop("The argument supplied for remote.jags must be a character string specifying the path to JAGS on the snow nodes, or a function (with no arguments) returning a character string specifying the path to JAGS (or not specified in which case JAGS will be found using findjags)", call.=FALSE)
		}
		if(strmethod=="rjparallel"){
			if(is.na(ncores)){
				ncores <- 2
				warning("Unable to detect the available number of cores on your machine - using a maximum of 2 cores as a default")
			}
			n.sims <- min(n.chains, ncores)
			if(!is.na(inputsims)){
				if(inputsims!=as.integer(inputsims) || inputsims <1) warning("The supplied n.sims option is not a positive integer and was ignored") else n.sims <- min(n.chains, inputsims)
			}
			progress.bar <- 'none'  # no point having a progress bar
			method <- runjags.rjparallel			
			extramethodargs <- c(extramethodargs,list(monitor=monitor, modules=modules, factories=factories, adapt=adapt, burnin=burnin, sample=sample, n.chains=n.chains, thin=thin, by=by, progress.bar=progress.bar))
			writefiles <- FALSE
		}
		if(strmethod=="background"){
			n.sims <- 1
			method <- runjags.background
		}
		if(strmethod=="bgparallel"){
			if(is.na(ncores)){
				ncores <- 2
				warning("Unable to detect the available number of cores on your machine - using a maximum of 2 cores as a default")
			}
			n.sims <- min(n.chains, ncores)
			if(!is.na(inputsims)){
				if(inputsims!=as.integer(inputsims) || inputsims <1) warning("The supplied n.sims option is not a positive integer and was ignored") else n.sims <- min(n.chains, inputsims)
			}
			method <- runjags.background
		}
		if(strmethod=='xgrid'){
			# Leave here as a failsafe in case extend.jags is used on an xgrid method object:
			test <- setup.xgrid(mgridpath=xgrid.options$mgridpath, hostname=xgrid.options$hostname, password=xgrid.options$password, testonly=TRUE)
			n.sims <- min(n.chains, xgrid.options$max.threads)
			method <- runjags.xgrid
		}
		if(strmethod=='rjags'){
			n.sims <- 1
			# n.sims can never be more than 1 for rjags in case we don't specify parallel methods!!!
			method <- runjags.rjags
			extramethodargs <- c(extramethodargs,list(monitor=monitor, modules=modules, factories=factories, adapt=adapt, burnin=burnin, sample=sample, n.chains=n.chains, thin=thin, by=by, progress.bar=progress.bar))
			writefiles <- FALSE
		}
		
		method.options <- c(method.options, list(os=jags.status$os, libpaths=jags.status$libpaths, n.sims=n.sims, jobname=jobname, cl=cl, remote.jags=remote.jags, rjags=c(list(rjags=rjags),extramethodargs)))
		method.options$jags <- jags.status$JAGS.path
		
		if(strmethod=='xgrid'){
			method.options <- c(method.options, xgrid.options[names(xgrid.options)!="jobname"])
			if(!identical(names(method.options), c("jags", "silent.jags", "jags.refresh", "batch.jags", "os", "libpaths", "n.sims", "jobname", "cl", "remote.jags", "rjags", "command", "customart", "jagspath", "submitandstop", "max.threads", "mgridpath", "hostname", "password"))) stop("Invalid method.options list provided - ensure xgrid jags jobs are started with the correct xgrid.xx.jags functions", call.=FALSE)
		}else{
			stopifnot(all(c("jags", "silent.jags", "jags.refresh", "batch.jags", "os", "libpaths", "n.sims", "jobname", "cl", "remote.jags", "rjags") %in% names(method.options)))	
		}
		
		stopifnot(is.function(method))
		
	}else{
		
		if(!class(method.options)=="list" || any(names(method.options)=="")) stop("The method.options argument must be provided as a named list")
		
		if(!any(names(method.options)=="n.sims")){
			warning("No 'n.sims' value provided in the method.options list; assuming a single simulation directory is required")
			method.options$n.sims <- 1
		}else{
			if(method.options$n.sims!=as.integer(method.options$n.sims) || method.options$n.sims <1){
				warning("The supplied n.sims option is not a positive integer and was ignored; assuming a single simulation directory is required")
				method.options$n.sims <- 1
			}else{
				method.options$n.sims <- min(n.chains, method.options$n.sims)
			}
		}
		
		if(tempdir){
			temp.directory <- tempfile('runjagsdir')
			jobname <-  temp.directory
			dir.create(temp.directory)
		}else{
			# Change jobname to match folder name:
			jobname <- new_unique('runjagsfiles', touch=TRUE, type='folder')
			if(jobname=="Directory not writable"){
				stop("Directory not writable", call.=FALSE)
			}
			temp.directory <- file.path(getwd(), jobname)
		}
		if(any(names(method.options)=="write.files")) write.files <- method.options$write.files
			
		strmethod <- "custom"
		
	}	

	if(method.options$n.sims != as.integer(method.options$n.sims)) stop("The number of simulations required must be an integer")
	n.sims <- method.options$n.sims
	
	# Check that the function matches the function arguments:
	if(any(is.na(pmatch(names(method.options), names(formals(method))))))
		warning("One or more method.options supplied does not match an argument to the supplied method function - these arguments will be ignored")
  
  stopifnot(length(method.options)==length(unique(names(method.options))))

	method.options <- method.options[!is.na(pmatch(names(method.options), names(formals(method))))]
		
	
	if(any(c("pd","full.pd","popt","dic","ped") %in% monitor)){
		if(n.chains < 2) stop("The DIC, PED, pD, full.pD and pOpt cannot be assessed with only 1 chain")
		if(strmethod%in%runjagsprivate$parallelmethod || n.sims > 1) stop("The DIC, PED, pD, full.pD and pOpt cannot be assessed when using parallel or separate chains")
	}

	sim.chains <- matrix(NA, ncol=ceiling(n.chains/n.sims), nrow=n.sims)
	sim.chains[1:n.chains] <- 1:n.chains

	sim.chains <- lapply(1:n.sims, function(x) return(sim.chains[x,!is.na(sim.chains[x,])]))
	nsim.chains <- sapply(sim.chains, length)

	stopifnot(sum(nsim.chains)==n.chains)
	
	method.options$rjags$sim.chains <- sim.chains
	
	# Only required for rjparallel:
	if(strmethod%in%runjagsprivate$rjagsmethod){
		method.options$rjags$origdatanames <- names(list.format(data))
	}
	
	real.runs <- as.integer(updates)
	ini.runs <- as.integer(burnin)
	adapt.runs <- as.integer(adapt)
	
	realmonitor <- monitor[! (monitor %in% c("dic","ped","popt","pd","full.pd")) ]
	# We have already checked for monitors but that included dic/popt/pd/full.pd - but we can't have zero real monitors (not including popt/full.pd/pd):
	if(length(realmonitor)==0) stop("Cannot run a model with only popt, pd or full.pd monitored - add a named variable (e.g. 'deviance') to the monitors", call.=FALSE)
	monitorcollapse <- paste(", thin(", thin, ")\nmonitor ", sep="")
	monitors <- paste("monitor ", paste(realmonitor, collapse=monitorcollapse), ", thin(", thin, ")\n", sep="")

	n.params <- length(monitor)
	params.names <- monitor
	
	total.updates <- n.params * updates
	
	results <- list()
	
	save.directory <- getwd()


	on.exit(setwd(save.directory))

	
  if(length(modules)==0)
    modules <- ''
  
	initstring <- paste(inits, "\n", sep="")
	
	resetbinpath=resetsyspath <- FALSE
	if(writefiles){
		
		setwd(temp.directory)

		cat(model, file="model.txt",sep="")
		close(file("model.txt"))
		cat(data, file="data.txt",sep="")  
		close(file("data.txt"))
		
		save(sim.chains, file="simchainsinfo.Rsave")
		close(file("simchainsinfo.Rsave"))
			
		for(s in 1:n.sims){

			sim <- paste("sim", s, sep=".")
			#system(paste('mkdir ', sim, sep=''))
			dir.create(sim)

			scriptstring <- ""
			
			if(!identical(modules,"")) for(i in 1:length(modules)){
				if(modules[[i]][1]=="runjags") stop("The runjags module is only available using the rjags method; to use the functions provided with other methods install (and specify using the module argument) the 'paretoprior' standalone module")
				scriptstring <- paste(scriptstring, if(modules[[i]][2]=='FALSE') "un", "load ", modules[[i]][1], "\n", sep="")
			}
			
			if(!identical(factories,"")) for(i in 1:length(factories)){
				scriptstring <- paste(scriptstring, "set factory \"", factories[[i]][1], "\" ", if(factories[[i]][3]=='TRUE') "on" else "off", ", type(", factories[[i]][2], ")\n", sep="")
			}

			scriptstring <- paste(scriptstring, "model in \"model.txt\"\n", sep="")
			if(data!=""){
				scriptstring <- paste(scriptstring, "data in \"data.txt\"\n", sep="")
			}

			scriptstring <- paste(scriptstring, "compile, nchains(", as.integer(nsim.chains[s]), ")\n", sep="")
			for(c in 1:nsim.chains[s]){
				i <- sim.chains[[s]][c]
				if(!is.na(inits[i]) && inits[i]!="") scriptstring <- paste(scriptstring, "parameters in \"inits", i, ".txt\", chain(", c, ")\n", sep="")
			}

			scriptstring <- paste(scriptstring, "initialize\n", sep="")
			
			####  TODO adapt and adapt.incomplete option needs sorting for JAGS 3 vs 4 and batch.jags FALSE/TRUE
			# To check adaptation I either do: 
			# adapt 0 -> error if it needs adapting, not otherwise
			# update 0 -> warning if adaptation incomplete
			# adapt 1 will cause an error in JAGS 3 but not 4 - needs sorting out properly		
			if(adapt.runs > 0){
				scriptstring <- paste(scriptstring, "adapt ", adapt.runs, "\n", sep="")
			}else{
				# If we have specified a burnin but adapt==0, allow the adapt to happen during burnin
				if(runjags.getOption('adapt.incomplete')=='error' && ini.runs==0){
					scriptstring <- paste(scriptstring, "adapt 0\n", sep="")
				}else{
					if(ini.runs == 0)   # Don't want 2 lots of update 0!
						scriptstring <- paste(scriptstring, "update ", ini.runs, "\n", sep="")
				}
			}			
			if(ini.runs > 0)
				scriptstring <- paste(scriptstring, "update ", ini.runs, "\n", sep="")
			####
			
			# trace monitor for deviance is just a standard monitor name (and not used to calculate DIC)
			scriptstring <- paste(scriptstring, if(any(c("dic","ped")%in%monitor)) paste("monitor deviance, type(mean) thin(", thin, ")\n", "monitor pD, type(mean) thin(", thin, ")\n", "monitor popt, type(mean) thin(", thin, ")\n", sep=""), monitors, if(any(monitor=="full.pd")) paste("monitor pD, thin(", thin, ")\n", sep=""), sep="")
			
			if(real.runs > 0) scriptstring <- paste(scriptstring, "update ", real.runs, "\n", sep="")
		
			for(c in 1:nsim.chains[s]){
				i <- sim.chains[[s]][c]
				scriptstring <- paste(scriptstring, "parameters to \"out", i, ".Rdump\", chain(", c, ")\n", sep="")
			}

			# Write out the mean pd, popt and deviances and then delete them so they dont get re-written by coda *
			if(any(c("dic","ped")%in%monitor))
				scriptstring <- paste(scriptstring, "coda deviance, stem(\"sim.", s, "/deviance\")\n", "coda pD, stem(\"sim.", s, "/pd\")\n", "coda popt, stem(\"sim.", s, "/popt\")\n", "monitor clear deviance, type(mean)\nmonitor clear pD, type(mean)\nmonitor clear popt, type(mean)\n", sep="")
			
			scriptstring <- paste(scriptstring, "coda *, stem(sim.", s, "/CODA)\n", sep="")
			
			scriptstring <- paste(scriptstring, "samplers to sim.", s, "/samplers.csv\n", sep="")	
			
			# model clear is used to detect the model being finished, update 0 is used to detect the model not having crashed
			scriptstring <- paste(scriptstring, "update 0\nmodel clear\nexit\n", sep="")

			output <- file(paste("sim.",s,"/script.cmd", sep=""), 'w')
			cat(scriptstring, file=output,sep="")  
			close(output)
			
		}

		for(i in 1:chains){
			cat(initstring[i], file=paste("inits", i, ".txt", sep=""),sep="")
			close(file(paste("inits", i, ".txt", sep="")))
		}


		os <- .Platform$OS.type	
		if(os == "windows"){		
			currentsyspath <- Sys.getenv('PATH')
			if(!grepl(method.options$libpaths$PATH,currentsyspath,fixed=TRUE)){
				Sys.setenv(PATH=paste(currentsyspath, ';', method.options$libpaths$PATH, sep=''))
				resetsyspath <- TRUE
			}

			currentsysbinpath <- Sys.getenv('LTDL_LIBRARY_PATH')
			if(!grepl(method.options$libpaths$LTDL_LIBRARY_PATH,currentsysbinpath,fixed=TRUE)){
				Sys.setenv(LTDL_LIBRARY_PATH=paste(currentsysbinpath, if(currentsysbinpath!='') ';', method.options$libpaths$LTDL_LIBRARY_PATH, sep=''))
				resetbinpath <- TRUE
			}		
		}	
	}
		
	result <- try(do.call(method, method.options), silent=TRUE)
	
	if(resetsyspath) Sys.setenv(PATH=currentsyspath)
	if(resetbinpath) Sys.setenv(LTDL_LIBRARY_PATH=currentsysbinpath)
	
#	fp <- c('sim.1/CODAchain1.txt', 'sim.1/CODAindex.txt')
#	repeat{
#		print(isOpen(textConnection(fp[1])))
#		print(isOpen(textConnection(fp[2])))
#		fi <- file.info(fp)
#		print(fi)
#		Sys.sleep(0.5)
#	}
	
	if(class(result)=="try-error"){
		failedjagsmodel <- model
		class(failedjagsmodel) <- "runjagsmodel"
		assign("model", failedjagsmodel, envir=failedjags)

		failedd <- data
		class(failedd) <- "runjagsdata"
		assign("data", failedd, envir=failedjags)
		
		failedi <- inits
		class(failedi) <- "runjagsinits"
		assign("inits", failedi, envir=failedjags)
		
		failedo <- as.character(result)
		class(failedo) <- "rjagsoutput"
		assign("output",failedo, envir=failedjags)
		
		stop(paste("The following error was encountered while attempting to run the JAGS model:  \n   ", gsub('Error : ','',as.character(result),fixed=TRUE), sep=""),call.=FALSE)
	}
		
	if(class(result)=="list"){
		results <- result
		if(any(names(results)=="complete")){
			result <- results$complete
			if(!all(names(results)==c('complete','mcmc','deviance.table','deviance.sum','pd','end.state','samplers')))
				stop('Unsupported return type from the JAGS method function - it should be a list with the elements complete, mcmc, deviance.table (may be NA), deviance.sum (may be NA), pd (may be NA), end.state, and samplers',call.=FALSE)
		}else{
			result <- results[[1]]
			names(results)[1] <- "complete"
		}
	}else{
		results <- list(complete=result)
	}
	
	if(class(result)=="character"){
		
		failedjagsmodel <- model
		class(failedjagsmodel) <- "runjagsmodel"
		assign("model", failedjagsmodel, envir=failedjags)

		failedd <- data
		class(failedd) <- "runjagsdata"
		assign("data", failedd, envir=failedjags)
		
		failedi <- inits
		class(failedi) <- "runjagsinits"
		assign("inits", failedi, envir=failedjags)
		
		failedo <- result
		class(failedo) <- "rjagsoutput"
		assign("output",failedo, envir=failedjags)

		stop(result)
	}
		
	if(!is.logical(result)) stop("Unsupported return type from the JAGS method function",call.=FALSE)
	
		# Done by on.exit call earlier:
	# setwd(save.directory)	
	
	results <- c(results, list(jobname=jobname, directory=normalizePath(temp.directory, winslash='/'), n.sims=n.sims, startedon=Sys.time()))
	return(results)		
	
}

				