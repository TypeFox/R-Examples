xgrid.start <- function(method, f, niters, object.list, file.list, threads, Rbuild, arguments, tempdir, keep.files, show.output, ...){
	
	if(class(object.list)=='character'){
		new.list <- vector('list')
		for(i in object.list){
			new.list <- c(new.list, list(get(i, envir=parent.frame(2))))
		}
		names(new.list) <- object.list
		object.list <- new.list
	}
	
	if(class(object.list)!='list') stop('The object.list supplied must be a named list or character vector', call.=FALSE)
	if(is.null(names(object.list)) & length(object.list) > 0) stop('The object.list supplied must be a named list or character vector', call.=FALSE)
	
	objects <- object.list
	
	xg.extra.args <- list(...)
	xg.arguments <- arguments
	
	if(any(names(objects)=="")) stop('All objects to be passed to the function (either supplied in the object.list or specified in the xgrid.run function call) must be named', call.=FALSE)
	if(any(unlist(lapply(names(objects), function(x) return(strsplit(x, 'xg.')[[1]][1])))=="")) stop("Sorry - you can't use any variable names with the prefix 'xg.' in the object list as they could cause naming conflicts with variables used to start and run the function on Xgrid", call.=FALSE)
		
	if(threads > niters){
		warning('The number of threads was reduced to match the number of iterations', call.=FALSE)
		threads <- niters
	}
			
	if(threads==1){
		method$separate.tasks <- FALSE
		method$separate.jobs <- FALSE
	}
	separatetasks <- method$separate.tasks
	separatejobs <- method$separate.jobs	
	
	if(class(objects)!='list') stop('The objects must be supplied as a named list', call.=FALSE)
	if(is.null(names(objects)) & length(objects)!=0) stop('The objects must be supplied as a named list', call.=FALSE)
	if(any(names(objects)=='')) stop('The objects must be supplied as a named list', call.=FALSE)
	
	if(length(objects)==0) objects$blank <- 0
	
	# Do this in a sub-function to avoid naming conflicts:
	saveobjects <- function(objects){
		for(i in 1:length(objects)){
			assign(names(objects)[i], objects[[i]])
		}
		save(list=names(objects), file='Robjects.Rsave')
		return(file.exists('Robjects.Rsave'))
	}
		
	gottoend <- FALSE
	interrupt <- FALSE
	
	tryCatch({
	
		save.directory <- getwd()
		on.exit(setwd(save.directory))

		if(tempdir){
			temp.directory <- tempfile('xgridrundir')
			dir.create(temp.directory)
		}else{
			method$jobname <- new_unique(method$jobname, touch=TRUE, type='folder')
			if((method$jobname=="Directory not writable")==TRUE){
				stop("Directory not writable", call.=FALSE)
			}
			temp.directory <- paste(getwd(), '/', method$jobname, sep='')
		}
		
		if(length(file.list)>0){
			for(p in file.list){
				success <- file.copy(from=p, to=temp.directory)
				if(!success) stop(paste("File '", p, "' not found in the working directory", sep=""), call.=FALSE)
			}
		}
		
		setwd(temp.directory)

		success <- saveobjects(objects)
		
		filesize <- as.numeric(system("ls -l Robjects.Rsave | awk '{print $5}'", intern=TRUE))
		if(filesize>method$max.filesize) stop('The object list supplied is greater than the maximum filesize specified.  You could reduce the number of objects supplied by making some large files available on a shared volume, or increase max.filesize (this might crash your xgrid controller...)', call.=FALSE)
		
		if(threads==niters){
			iterations <- as.list(1:niters)
		}else{
			t <- 1
			iterations <- vector('list', length=threads)
			for(i in 1:niters){
				iterations[[t]] <- c(iterations[[t]], i)
				t <- t+1
				if(t>threads) t <- 1
			}
		}
		
		# The maximum filesize is TOTAL, so assume all threads are similar and:
		if(separatetasks) max.filesize <- method$max.filesize/threads else max.filesize <- method$max.filesize
		
		xg.f <- f
		xg.iterations <- iterations
		xg.max.filesize <- max.filesize
		
		save(xg.f, xg.iterations, xg.arguments, xg.extra.args, xg.max.filesize, file='Rfunction.Rsave')
		
		cat(method$customart, "\n", file="customart.sh")
				
		# Validity of Rbuild is checked by setup.xgrid		
		if(is.na(Rbuild)) Rbuild <- ''
		if(Rbuild==' ') Rbuild <- ''				
		
		loadpackagestring <- paste("packagesuccess <- try({ ", if(length(method$packages)>0) paste("library(", paste(names(method$packages), collapse="); library("), "); ", sep=""), " }); ", sep="")
				
		if(separatejobs){
			jobnames <- paste(method$jobname, 1:threads, sep='.thread.')
			for(s in 1:threads){
				cat('#!/bin/sh
pid=$$

alias myr=\'', method$Rpath, '\'
', if(Rbuild!="") paste('if [ -f \'', method$Rpath, Rbuild, '\' ]; then
alias myr=\'', method$Rpath, Rbuild, '\'
fi
test=`\'', method$Rpath, Rbuild, '\' --version 2>&1`
if [ `echo $?` != 0 ]; then
alias myr=\'', method$Rpath, '\'
fi
', sep=''), '

( ( (myr --slave --no-save -e "xg.task <- ', s, '; xg.path <- Sys.getenv(\'PATH\'); xg.path <- paste(xg.path, \':/bin:/usr/local/bin:/usr/bin\', sep=\'\'); Sys.setenv(PATH=xg.path); load(\'Robjects.Rsave\', envir=.GlobalEnv); load(\'Rfunction.Rsave\', envir=.GlobalEnv); xg.iteration <- xg.iterations[[xg.task]]; ', loadpackagestring, 'for(i in 1:length(xg.iteration)){ xg.i <- xg.iteration[i]; if(class(packagesuccess)==\'try-error\'){ xg.temp <- \'Error loading R packages\' } else { xg.temp <- \'Error\'; try(xg.temp <- do.call(xg.f, c(list(xg.arguments[[xg.i]]), xg.extra.args)))}; assign(paste(\'iteration.\', xg.i, sep=\'\'), xg.temp); cat(\'<xgrid>{control = statusUpdate; percentDone = \', round(i/length(xg.iteration)*100, digits=0), \'; }</xgrid>\', sep=\'\')}; save(list=paste(\'iteration.\', xg.iteration, sep=\'\'), file=paste(\'results.\', xg.task, \'.Rsave\', sep=\'\')); xg.filesize <- as.numeric(system(paste(\'ls -l results.\', xg.task, \'.Rsave | awk \\\'{print \\$5}\\\'\', sep=\'\'), intern=TRUE)); if(xg.filesize>xg.max.filesize){ xg.temp <- \'Maximum file size exceeded.  You could either try again using more jobs, reduce the number of objects returned by the function, or save files to a shared drive rather than returning them\'; for(xg.i in xg.iteration) assign(paste(\'iteration.\', xg.i, sep=\'\'), xg.temp); save(list=paste(\'iteration.\', xg.iteration, sep=\'\'), file=paste(\'results.\', xg.task, \'.Rsave\', sep=\'\'))};"; echo $? > .retstat.$pid) 2>&1 1>&3 | tee Rerror.$1.txt) 3>&1 1>&2) 2>&1 | tee Rout.$1.txt

# This makes sure the process has finished before continuing:
wait

returnstat=`cat < .retstat.$$`
rm .retstat.$$

exit $returnstat

', sep='', file=jobnames[s])
				Sys.chmod(jobnames[s])
			}
		}else{
			cat('#!/bin/sh
pid=$$

alias myr=\'', method$Rpath, '\'
', if(Rbuild!="") paste('if [ -f \'', method$Rpath, Rbuild, '\' ]; then
alias myr=\'', method$Rpath, Rbuild, '\'
fi
test=`\'', method$Rpath, Rbuild, '\' --version 2>&1`
if [ `echo $?` != 0 ]; then
alias myr=\'', method$Rpath, '\'
fi
', sep=''), '

echo "" > Rout.', if(separatetasks) '$1' else '1', '.txt
echo "" > Rerror.', if(separatetasks) '$1' else '1', '.txt
', if(separatetasks) '( ( echo "\nTask "$1":" 2>&1 1>&3 | tee -a Rerror.$1.txt) 3>&1 1>&2) 2>&1 | tee -a Rout.$1.txt', '

( ( (myr --slave --no-save -e "xg.task <- ', if(separatetasks) '$1' else '1', '; xg.path <- Sys.getenv(\'PATH\'); xg.path <- paste(xg.path, \':/bin:/usr/local/bin:/usr/bin\', sep=\'\'); Sys.setenv(PATH=xg.path); load(\'Robjects.Rsave\', envir=.GlobalEnv); load(\'Rfunction.Rsave\', envir=.GlobalEnv); xg.iteration <- xg.iterations[[xg.task]]; ', loadpackagestring, 'for(i in 1:length(xg.iteration)){ xg.i <- xg.iteration[i]; if(class(packagesuccess)==\'try-error\'){ xg.temp <- \'Error loading R packages\' } else { xg.temp <- \'Error\'; try(xg.temp <- do.call(xg.f, c(list(xg.arguments[[xg.i]]), xg.extra.args)))}; assign(paste(\'iteration.\', xg.i, sep=\'\'), xg.temp); cat(\'<xgrid>{control = statusUpdate; percentDone = \', round(i/length(xg.iteration)*100, digits=0), \'; }</xgrid>\', sep=\'\')}; save(list=paste(\'iteration.\', xg.iteration, sep=\'\'), file=paste(\'results.\', xg.task, \'.Rsave\', sep=\'\')); xg.filesize <- as.numeric(system(paste(\'ls -l results.\', xg.task, \'.Rsave | awk \\\'{print \\$5}\\\'\', sep=\'\'), intern=TRUE)); if(xg.filesize>xg.max.filesize){ xg.temp <- \'Maximum file size per task exceeded.  You could either try again using multiple jobs, reduce the number of objects returned by the function, or save files to a shared drive rather than returning them\'; for(xg.i in xg.iteration) assign(paste(\'iteration.\', xg.i, sep=\'\'), xg.temp); save(list=paste(\'iteration.\', xg.iteration, sep=\'\'), file=paste(\'results.\', xg.task, \'.Rsave\', sep=\'\'))};"; echo $? > .retstat.$pid) 2>&1 1>&3 | tee -a Rerror.', if(separatetasks) '$1' else '1', '.txt) 3>&1 1>&2) 2>&1 | tee -a Rout.', if(separatetasks) '$1' else '1', '.txt

', if(separatetasks) '( ( echo "\n" 2>&1 1>&3 | tee -a Rerror.$1.txt) 3>&1 1>&2) 2>&1 | tee -a Rout.$1.txt', '

# This makes sure the process has finished before continuing:
wait

returnstat=`cat < .retstat.$$`
rm .retstat.$$

exit $returnstat

', sep='', file=method$jobname)
			Sys.chmod(method$jobname)
	
		}

		# Separatejobs will only ever be for xgrid.submit:
		if(separatejobs){

			jobnum <- integer(threads)

			for(s in 1:threads){

				cat('#!/bin/sh
cmd="', jobnames[s], '"
ntasks=1
job=', s, '
indir="', temp.directory, '"
', method$command, if(!show.output) ' >/dev/null', ' &
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
				if(file.exists("mgridpid.txt")){
					interrupt <- TRUE
					pid <- as.numeric(readLines('mgridpid.txt'))
					success <- system(paste("kill ", pid, sep=""))
					stop("The process was interrupted while submitting the xgrid job", call.=FALSE)
				}
		
				if(success!=0){
					stop("An unknown error occured while starting the xgrid command", call.=FALSE)
				}
				
				Sys.sleep(1)
				
				if(file.exists('jobid.txt')) tjobnum <- paste(readLines('jobid.txt'), collapse='\n') else tjobnum <- paste(readLines('.starterout.txt'), collapse='\n')

				jobnum[s] <- gsub('[^[:digit:]]', '', paste(tjobnum, collapse=''))

				if(jobnum[s]=='' | as.numeric(jobnum[s])>10^6) stop(paste("There was an error submitting job number ", s, " to Xgrid - no Job ID was returned by mgrid/xgrid", sep=''))
				
				swcat("Job ", s, " of ", threads, " submitted to xgrid\n", sep="")
				
				# Separatejobs will only ever be for xgrid.submit!
				
			}

		}else{

			cat('#!/bin/sh

cmd="', method$jobname, '"
ntasks=', threads, '
job=1
indir="', temp.directory, '"
', method$command, if(!show.output) ' >/dev/null', ' &
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
			if(file.exists("mgridpid.txt")){
				interrupt <- TRUE
				pid <- as.numeric(readLines('mgridpid.txt'))
				success <- system(paste("kill ", pid, sep=""))
				stop("The process was interrupted while submitting or waiting for the results of the xgrid job", call.=FALSE)
			}
		
			if(success!=0){
				stop("An unknown error occured while starting the xgrid command", call.=FALSE)
			}
				
			if(method$submitandstop){
				Sys.sleep(1)
				
				if(file.exists('jobid.txt')) jobnum <- paste(readLines('jobid.txt'), collapse='\n') else jobnum <- paste(readLines('starteroutput.txt'), collapse='\n')

				jobnum <- gsub('[^[:digit:]]', '', paste(jobnum, collapse=''))
			
				if(jobnum=='' | as.numeric(jobnum)>10^6) stop("There was an error submitting your job to Xgrid - no Job ID was returned")
			
				swcat('Your job (name: "', method$jobname, '";  ID: "', jobnum, '") has been succesfully uploaded to xgrid\n', sep='')

			}else{
				jobnum <- NA
				swcat("Finished running the job\n")
			}

		}
	
	
		if(method$submitandstop){
			cat(jobnum, file='jobid.txt', sep='\n')
			savelist <- ls()
			savelist <- savelist[savelist!='keep.files' & savelist!='cleanup' & savelist!='show.output' & savelist!='partialretrieve']
			save(list=savelist, file='workingobj.Rsave')
		}
	
		gottoend <- TRUE		
			
	}, finally={
		
		jagsrun <- FALSE
		directory <- temp.directory
		jobname <- method$jobname
		jobid <- jobnum
		save(jagsrun, directory, jobname, jobid, threads, niters, iterations, file="info.Rsave")
		
		if(!method$submitandstop & tempdir & keep.files){
			setwd(save.directory)
			new.directory <- new_unique(method$jobname, touch=TRUE, type='folder')
			if((new.directory=="Directory not writable")==TRUE){
				warning("Xgrid files could not be copied to the working directory as it is not writable")
			}else{
				file.copy(from=paste(temp.directory, list.files(temp.directory), sep=.Platform$file.sep), to=new.directory, recursive=TRUE)
				swcat("Xgrid files were saved to the '", new.directory, "' folder in your current working directory\n", sep="")
			}
			unlink(temp.directory, recursive=TRUE)
		}
		
		if(!gottoend){
			if(interrupt){
				stop("The process was terminated by the user", call.=FALSE)
			}else{
				stop("An unknown error occured while setting up or running the model", call.=FALSE)
			}
		}
	})
	
	return(list(directory=temp.directory, jobname=method$jobname, jobid=jobnum, threads=threads, niters=niters, iterations=iterations))
}

xgrid.read <- function(threads, niters, iterations){
	
	success <- try(results <- vector('list', length=length(unlist(iterations))), silent=TRUE)
	if(inherits(success, 'try-error')) stop('The results were not in the expected format - if this job was started using xgrid.submit.jags then use the xgrid.results.jags function to retrieve the results')
	names(results) <- paste('iteration.', 1:length(results), sep='')
	failed <- logical(length(unlist(iterations)))
	for(t in 1:length(iterations)){
		
		options(show.error.messages = FALSE)
		try({
		suppressWarnings(success <- try(load(paste('results.', t, '.Rsave', sep='')), silent=TRUE))
		})
		options(show.error.messages = TRUE)
		
		if(inherits(success, 'try-error')){
			for(i in iterations[[t]]){
				results[[paste('iteration.', i, sep='')]] <- 'Job or task incomplete'
				failed[i] <- TRUE
			}
		}else{
			for(i in iterations[[t]]){
				results[[paste('iteration.', i, sep='')]] <- get(paste('iteration.', i, sep=''))
			}
		}
	}
	
	#if(any(failed)) cat('Results for iterations ', paste(which(failed),collapse=','), ' were not returned\n', sep='')
#	if(any(failed)) swcat('\nJob ', round(1-((sum(failed)/length(failed))), digits=2)*100, '% complete\n', sep='') else swcat('\nJob complete\n')	
	# Don't think tacking this onto the end of the return is really appropriate...
	#results <- c(results, list(output=paste(xgridoutput,collapse='\n')))
	return(list(results=results, failed=failed))
	
	
}

xgrid.run <- function(f=function(iteration){}, niters=1, object.list=list(), file.list=character(0), max.threads=100, arguments=as.list(1:niters), Rversion="", packages=list(), artfun=function() writeLines("1"), email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, cleanup=TRUE, showprofiles=FALSE, Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), tempdir=FALSE, keep.files=FALSE, show.output=TRUE, threads=min(niters, max.threads), ...){
	
	method <- setup.xgrid(separate=FALSE, JAGSversion="", Rversion=Rversion, packages=packages, artfun=artfun, email=email, profiling=profiling, cpuarch=cpuarch, minosversion=minosversion, queueforserver=queueforserver, hostnode=hostnode, forcehost=forcehost, ramrequired=ramrequired, jobname=jobname, cleanup=cleanup, showprofiles=showprofiles, jagspath="", Rpath=Rpath, Rbuild=Rbuild, max.filesize=max.filesize, mgridpath=mgridpath, hostname=hostname, password=password, submitandstop=FALSE, jagsrun=FALSE)

	submission <- xgrid.start(method=method, f=f, niters=niters, object.list=object.list, file.list=file.list, threads=threads, Rbuild=Rbuild, arguments=arguments, tempdir=tempdir, keep.files=keep.files, show.output=show.output, ...)
	
	save.directory <- getwd()
	on.exit(setwd(save.directory))

	setwd(submission$directory)
	
	results <- xgrid.read(threads=submission$threads, niters=submission$niters, iterations=submission$iterations)
	
	if(any(results$failed)) warning("There was an unexpected error running the job(s);  one or more of the tasks failed to return a result", call.=FALSE)

	return(results$results)
}

xgrid.submit <- function(f=function(iteration){}, niters=1, object.list=list(), file.list=character(0), max.threads=100, arguments=as.list(1:niters), Rversion="", packages=list(), artfun=function() writeLines("1"), email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), show.output=TRUE, separate.jobs=FALSE, threads=min(niters, max.threads), ...){
	
	method <- setup.xgrid(separate=separate.jobs, JAGSversion="", Rversion=Rversion, packages=packages, artfun=artfun, email=email, profiling=profiling, cpuarch=cpuarch, minosversion=minosversion, queueforserver=queueforserver, hostnode=hostnode, forcehost=forcehost, ramrequired=ramrequired, jobname=jobname, cleanup=FALSE, showprofiles=FALSE, jagspath="", Rpath=Rpath, Rbuild=Rbuild, max.filesize=max.filesize, mgridpath=mgridpath, hostname=hostname, password=password, submitandstop=TRUE, jagsrun=FALSE)

	submission <- xgrid.start(method=method, f=f, niters=niters, object.list=object.list, file.list=file.list, threads=threads, Rbuild=Rbuild, arguments=arguments, tempdir=FALSE, keep.files=TRUE, show.output=show.output, ...)
	
	return(list(jobname=submission$jobname, directory=submission$directory, jobid=submission$jobid))
	
}


xgrid.results <- function(jobinfo, wait=TRUE, partial.retrieve=!wait, cleanup=!partial.retrieve, show.output=TRUE){
	
	# Special clause to see if jobname is actually a jags job (if jobname/jagsinfo.Rsave exists and jobname/jobname text file contains jagserror.txt and jagsout.txt) - bounce to xgrid.results.jags (with warning?):
	if(class(jobinfo)=="runjagsbginfo"){
		warning("Attempting to retrieve results from an Xgrid JAGS run using the non-JAGS Xgrid function - redirecting the function call to xgrid.results.jags...")
		return(xgrid.results.jags(background.runjags.object=jobinfo, wait=wait, cleanup=cleanup))			
	}
	
	if(wait & partial.retrieve) stop("Can't use partial.retrieve option if waiting for the xgrid job to finish!", call.=FALSE)
	
	swcat("Retrieving xgrid results...\n\n")
	
	# First go to xgrid.retrieve which checks the jobname exists and gets files back from xgrid (silent is always FALSE as if silent.jags was specified, then jags output is binned):
	output <- xgrid.retrieve(jobinfo=jobinfo, wait=wait, silent=!show.output, cleanup=cleanup, partialretrieve=partial.retrieve, jags=FALSE)
	
	save.directory <- getwd()
	on.exit(setwd(save.directory))

	setwd(output$directory)
	
	# To prevent warnings about no visible bindings:
	jagsrun = directory = jobname = jobid = threads = niters = iterations <- NA
	# In xgrid.submit and xgrid.jags.submit (minus threads,niters,iterations):
	# save(jagsrun, directory, jobname, jobid, threads, niters, iterations, file="info.Rsave")	
	load('info.Rsave')	
	
	results <- xgrid.read(threads=threads, niters=niters, iterations=iterations)

	if(any(results$failed)) swcat('Results for iterations ', paste(which(results$failed),collapse=','), ' were not returned\n', sep='')
	if(any(results$failed)) swcat('\nJob ', round(1-((sum(results$failed)/length(results$failed))), digits=2)*100, '% complete\n', sep='') else swcat('\nJob complete\n')	

	return(results$results)

}

xapply <- function(X, FUN, method.options=list(), ...){
	
	if(class(X)!='list') X <- as.list(X)
	arguments <- method.options
	arguments$f <- as.function(FUN)
	arguments$niters <- length(X)
	arguments$arguments <- X
	
	arguments <- c(arguments, list(...))
	
	return(do.call(xgrid.run, arguments, quote=FALSE))
	
}
