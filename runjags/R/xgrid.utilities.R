xgrid.jobs <- function(comment=FALSE, user=FALSE, jobs=10, mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"))
{
	
	# Call setup.xgrid just to make sure xgrid is available etc
	test <- setup.xgrid(mgridpath=mgridpath, hostname=hostname, password=password, testonly=TRUE)
	
	command <- paste(mgridpath, " -l ", if(comment) "-m ", if(user) "-u ", "-j ", jobs, " 2>&1", sep="")

	# Note ugly workarounds to stop console appearing to freeze on Aqua:
	if(.Platform$GUI=="AQUA"){
		swcat("Retrieving current job list from xgrid...")
		system("printf '\n'")
		wrapper <- tempfile()
		cat("#!/bin/bash\n", command, "\nexit 0\n", sep="", file=wrapper)
		Sys.chmod(wrapper)
		output <- system(wrapper, intern=TRUE)
		unlink(wrapper)
		swcat(gsub("Retrieving current job list from xgrid...\r","",paste(output,collapse="\n"),fixed=TRUE), "\n")
		
	}else{
		system(command, intern=FALSE)
	}
		
	invisible(command)
}

xgrid.delete <- function(jobinfo, keep.files=FALSE){
	
	swcat("Deleting xgrid job...\n")
	
	# Note we ignore jobinfo$jobid - jobnum is saved in the workingobj.Rsave file:
	if(class(jobinfo)=="list"){
		directory <- jobinfo$directory
		foldername <- jobinfo$jobname	
		jobnum <- jobinfo$jobid	
	}else{
		directory <- NA
		foldername <- jobinfo
		jobnum <- NA
	}
	if(class(foldername)!="character"){
		stop("The jobinfo must be supplied either as a character containing a path to the submitted folder, or the return value of an xgrid submit function", call.=FALSE)
	}
	
	# Find the relevant folder - first look for a matching folder in the working directory, and then go for the directory argument for an absolute path:
	if(!file.exists(foldername)){		
		if(!is.na(directory) && file.exists(directory)){
			foldername <- directory
		}else{
			stop("The folder specified by jobinfo does not exist", call.=FALSE)
		}
	}
	
	save.directory <- getwd()
	on.exit(setwd(save.directory))
	setwd(foldername)

	# Read the jobid and check it matches:
	if(!file.exists("info.Rsave")) stop("The required info.Rsave file was not found in the specified folder; ensure you have specified the path to a folder created by the xgrid submit functions", call.=FALSE)
	
	# To prevent warnings about no visible bindings:
	jagsrun = directory = jobname = jobid = threads = niters = iterations <- NA
	# In xgrid.submit and xgrid.jags.submit (minus threads,niters,iterations):
	# save(jagsrun, directory, jobname, jobid, threads, niters, iterations, file="info.Rsave")	
	load('info.Rsave')
	
	# Note we ignore jobinfo$jobid - jobnum is saved in the workingobj.Rsave file:
	if(!is.na(jobnum) && jobnum!=jobid) jobnum <- jobid
	
	if(length(jobnum) > 1) separatejobs <- TRUE else separatejobs <- FALSE
	n.sims <- length(jobnum)
		
	if(separatejobs){
		
		for(s in 1:n.sims){
			xgriddeleteout <- system(paste('xgrid -job delete -id ', jobnum[s], sep=''), intern=TRUE)
			if(paste(xgriddeleteout, collapse='')!='{}') warning(paste('Possible error deleting xgrid job number ', jobnum[s], ' - please check this manually', sep=''))
		}
		
	}else{
		
		xgriddeleteout <- system(paste('xgrid -job delete -id ', jobnum, sep=''), intern=TRUE)
		if(paste(xgriddeleteout, collapse='')!='{}') warning(paste('Possible error deleting xgrid job number ', jobnum, ' - please check this manually', sep=''))
	}
	
	setwd(save.directory)
	
	if(!keep.files){
		unlink(foldername, recursive=TRUE)
	}
	
	swcat("Job deleted\n")
	
	invisible(TRUE)
}


setup.xgrid <- function(separate=FALSE, JAGSversion=">=2.0.0", Rversion="", packages=list(), artfun=function() writeLines("1"), email=NA, profiling=TRUE, cpuarch=NA, minosversion=NA, queueforserver=FALSE, hostnode=NA, forcehost=FALSE, ramrequired=10, jobname=NA, cleanup=TRUE, showprofiles=FALSE, jagspath='/usr/local/bin/jags', Rpath='/usr/bin/R', Rbuild='64', max.filesize="1GB", mgridpath=system.file("xgrid", "mgrid.sh", package="runjags"), hostname=Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), password=Sys.getenv("XGRID_CONTROLLER_PASSWORD"), submitandstop=FALSE, jagsrun=TRUE, testonly=FALSE){
	
	if(showprofiles & !profiling){
		warning("Can't use showprofiles without profiling - setting showprofiles to FALSE")
		showprofiles <- FALSE
	}
	
	# Check xgrid is available
	if(.Platform$OS.type=='windows'){
		stop('Xgrid functions are only available on machines running Mac OS X 10.5 (Leopard) or later and with access to an Xgrid controller', call.=FALSE)
	}
	
	xgridavail <- suppressWarnings(system('xgrid 2>&1', intern=TRUE))
	if(length(xgridavail)==1){
		stop('Xgrid is not available on this machine.  Xgrid functions are only available on machines running Mac OS X 10.5 (Leopard), OS X 10.6 (Snow Leopard) or OS X 10.7 (Lion) and with access to an Xgrid controller', call.=FALSE)
	}
		
	osvers <- as.numeric(strsplit(system("sysctl -n kern.osrelease", intern=TRUE), ".", fixed=TRUE)[[1]][1])-4
	if(osvers < 5){
		stop('Xgrid is not available on this machine.  Xgrid functions are only available on machines running Mac OS X 10.5 (Leopard), OS X 10.6 (Snow Leopard) or OS X 10.7 (Lion) and with access to an Xgrid controller', call.=FALSE)
	}
	
	if(hostname!=Sys.getenv("XGRID_CONTROLLER_HOSTNAME")){
		Sys.setenv(XGRID_CONTROLLER_HOSTNAME=hostname)	
	}else{
		if(hostname=="") stop("Xgrid hostname not specified.  Note that system environmental variables set in the .profile file are not accessible from the GUI version of R - either specify the hostname and password to xgrid.run, export the XGRID_CONTROLLER_HOSTNAME and password variables using Sys.setenv(), or use the console version of R instead")
	}
	if(password!=Sys.getenv("XGRID_CONTROLLER_PASSWORD")){
		Sys.setenv(XGRID_CONTROLLER_PASSWORD=password)	
	}
	
	xgridoutput <- suppressWarnings(system('expect -c "# exp_internal 1" -c "spawn xgrid -job list" -c"expect {Password:  { send \r\n; interact }eof { exit }}" 2>&1', intern=TRUE))
	xgridoutput <- paste(gsub("\r","\n",xgridoutput),collapse="")
	
	if(grepl("Unable to connect",xgridoutput)){
		stop(paste("The Xgrid controller '", Sys.getenv("XGRID_CONTROLLER_HOSTNAME"), "' could not be found.  Make sure the environmental variables XGRID_CONTROLLER_HOSTNAME (and, if required, XGRID_CONTROLLER_PASSWORD) are set to the match the details for the controller you are trying to connect to, or specify the arguments to the xgrid functions.", sep=""), call.=FALSE)
	}
	if(grepl("Authentication failed",xgridoutput) & grepl("Password:",xgridoutput)){
		stop(paste("The blank Xgrid controller password specified is incorrect.  Make sure the environmental variable XGRID_CONTROLLER_PASSWORD is set to the match the details for the controller you are trying to connect to, or specify the argument to the xgrid functions.", sep=""), call.=FALSE)
	}
	if(grepl("Authentication failed",xgridoutput) & !grepl("Password:",xgridoutput)){
		stop(paste("The Xgrid controller password '", Sys.getenv("XGRID_CONTROLLER_PASSWORD"), "' is incorrect.  Make sure the environmental variable XGRID_CONTROLLER_PASSWORD is set to the match the details for the controller you are trying to connect to, or specify the argument to the xgrid functions.", sep=""), call.=FALSE)
	}
	if(!grepl("Authentication failed",xgridoutput) & grepl("Password:",xgridoutput)){
		stop("The password on the controller is set as blank and Xgrid is continually prompting for a password, which will cause problems with using this script.  You could try updating to the latest version of OS X, or just set a password on your controller.", call.=FALSE)
	}
	
	if(!file.exists(mgridpath)) stop("Supplied path to mgrid script is invalid", call.=FALSE)
	
	# Bail here if we're only testing Xgrid functions (xgrid.jobs):
	if(testonly) return(TRUE)
	
	if(!any(Rbuild==c("64", "32", ""))) warning(paste('Ignoring unrecognised Rbuild "', Rbuild, '".  Use one of "32", "64" or leave blank if no preference', sep=""), call.=FALSE)
		
	if(class(max.filesize)=="numeric" | class(max.filesize)=="integer"){
		max.filesize <- max.filesize * 1024^3# DEFAULT IS GB
	}else{
		if(class(max.filesize)!="character") stop("max.filesize must be either a numeric or character value")
		str.time <- strsplit(max.filesize, "")[[1]]

		time.unit <- suppressWarnings(str.time[is.na(as.numeric(str.time)!=str.time)])
		time.unit <- tolower(time.unit[time.unit!=" "][1])
		max.filesize <- suppressWarnings(as.numeric(paste(na.omit(str.time[as.numeric(str.time)==str.time]) ,collapse="")))
	
		max.filesize <- max.filesize * switch(time.unit, g=1024^3, m=1024^2, k=1024, b=1, NA)
		if(is.na(max.filesize)) stop("Unrecognised unit of max file size -'", time.unit, "'")
	}
	
	artfails <- FALSE
	if(class(artfun)!="function"){
		stop("Non function argument supplied to artfun", call.=FALSE)
	}else{
		if(!is.null(formals(artfun))){
			stop("The function supplied to artfun must take exactly 0 parameters!", call.=FALSE)
		}
		output <- suppressWarnings(as.numeric(gsub("[1] ","",capture.output(artfun()),fixed=TRUE)))
		if(length(output)!=1 || is.na(output)) stop("Supplied artfun function does not cat() a numeric value as expected (ensure required libraries are loaded on this machine)", call.=FALSE)
		if(output==0){
			warning("Supplied artfun function fails for this machine!", call.=FALSE)
			artfails <- TRUE
		}
	}
	
		
	# If jagsrun separate means separatetasks (separatejobs not available), if !jagsrun it means separatejobs rather than separate tasks (xgrid.start function sets both to FALSE if 1 thread):
	if(jagsrun){
		separate.tasks <- separate
		separate.jobs <- FALSE
	}else{
		separate.tasks <- !separate
		separate.jobs <- separate			
	}
	
	# Create mgrid command here:
	
	if(is.na(jobname)){
		if(jagsrun) jobname <- paste('xgridrunJAGS.', Sys.info()['user'], sep='') else jobname <- paste('xgridrunR.', Sys.info()['user'], sep='')		
	}
	
	if(forcehost & is.na(hostnode)) stop("You must specify one or more (separated by ':') hostnodes to be able to use the forcehost argument", call.=FALSE)
	if(!is.na(cpuarch)){
		cpuarch <- tolower(cpuarch)
		if(!any(tolower(cpuarch)==c("intel","i386","ppc","x86_64","intel32","intel64"))) stop("Unrecognised cpuarch type: specify one of 'intel' (will use 32 or 64 bit machines), 'i386' (32 bit intel only), 'x86_64' (64 bit intel only), or 'ppc' (will use 32 or 64 bit ppc machines)", call.=FALSE)
	}
	
	rankscript <- ""
	customrank <- FALSE
	if(class(profiling)=="character"){
		customrank <- TRUE
		rankscript <- profiling
		if(!file.exists(rankscript)) stop("Specified node profiling script does not exist", call.=FALSE)
		# Check and read file in later on:
		profiling <- TRUE
		if(jagsrun){
			warning("The existence of JAGS is not automatically checked when a custom profiling script is provided", call.=FALSE)
		}else{
			warning("The existence of R (and any specified R packages) is not automatically checked when a custom profiling script is provided", call.=FALSE)
		}
	}
#	if(!profiling){
#		if(jagsrun){
#			warning("The existence of JAGS is not automatically checked when profiling is disabled", call.=FALSE)
#		}else{
#			warning("The existence of R (and any specified R packages) is not automatically checked when profiling is disabled", call.=FALSE)
#		}
#	}
	
	# Cleanup and showprofiles not used if submitandstop - should never get this error when called from parent functions anyway:
	if(submitandstop & (cleanup | showprofiles)) stop("Can't use submitandstop and cleanup or showprofiles options together")
	
	if(is.na(hostnode) & forcehost) stop("You must specify one or more valid hostnode names when using the forcehost option")
	
	# Check minosversion is numeric:
	if(!is.na(minosversion)){
		suppressWarnings(minosversion <- as.numeric(gsub("[[:alpha:]]", "", minosversion)))
		if(is.na(minosversion) || length(minosversion)==0){
			stop("Unrecognised minimum Mac OS X version supplied", call.=FALSE)
		}
	}
	
	mgridcmd <- paste(system.file("xgrid", "mgrid.sh", package="runjags"), " -n ", jobname, " ", if(!submitandstop) "-w ", if(cleanup) "-x ", if(showprofiles) "-p ", if(!is.na(email)) paste("-e ", email, " ", sep=""), if(!is.na(hostnode)) paste("-h ", paste('"', hostnode, '"', collapse='":"', sep=''), " ", sep=""), if(forcehost) "-f ", if(queueforserver) "-q ", if(!is.na(cpuarch)) paste("-c ", cpuarch, " ", sep=""), if(!is.na(minosversion)) paste("-k ", minosversion, " ", sep=""), "-r ", ramrequired, " ", if(!profiling) "-z none ", if(customrank) "-z customart.sh " else "-a customart.sh ", '-i "$indir" -t $ntasks "$cmd"', sep="")
#	print(mgridcmd)
	
	artfile <- tempfile()
				
	if(jagsrun & !customrank){
		
		Sys.setenv(mgridcustomartname="JAGS_Score")
		
		# Make jagsart.sh file:
		if(is.na(JAGSversion)) JAGSversion <- ">=0"
		if(JAGSversion=="") JAGSversion <- ">=0"
	
		matchtype <- "-ge"
		JAGSversion <- gsub(">=","",JAGSversion)
		if(grepl(">",JAGSversion)) matchtype <- "-gt"
		if(grepl("<=",JAGSversion)) matchtype <- "-le"
		JAGSversion <- gsub("<=","",JAGSversion)
		if(grepl("<",JAGSversion)) matchtype <- "-lt"
		JAGSversion <- gsub("<","",JAGSversion)
		if(grepl("=",JAGSversion)) matchtype <- "-eq"
		JAGSversion <- gsub("=","",JAGSversion)
	
		num.version <- numeric(3)
		numout <- as.numeric(strsplit(JAGSversion, ".", fixed=TRUE)[[1]])
		num.version[1:length(numout)] <- numout
		num.version <- num.version * c(10^6,10^3,1)
		requiredvers <- num.version[1]+num.version[2]+num.version[3]
		
		art <- cat("#!/bin/bash
# For some reason which (and /usr/bin/which) seem to fail on private controllers at least, so use:
if [ ! -f '", jagspath, "' ]; then
	echo 0
else
	vers=`echo 'exit' | '", jagspath, "' | awk 'NR > 1 { exit }; 1' | awk '{print $4}' | awk -F . '{print $1*1000000 + $2*1000 + $3 }'`
	if [ $vers ", matchtype, " ", format(requiredvers, scientific=FALSE), " ]; then
		echo 1
	else
		echo 0
	fi
fi
exit 0
", sep="", file=artfile)
	
	}
	if(!jagsrun & !customrank){
		
		Sys.setenv(mgridcustomartname="R_Score")
		
		# Rversion - eg ">=2" - put note in help for Rbuild that you can ensure you get 64 vs 32 bit builds with an ART function
		if(is.na(Rversion)) Rversion <- ">=0"
		if(Rversion=="") Rversion <- ">=0"
	
		matchtype <- "-ge"
		Rversion <- gsub(">=","",Rversion)
		if(grepl(">",Rversion)) matchtype <- "-gt"
		if(grepl("<=",Rversion)) matchtype <- "-le"
		Rversion <- gsub("<=","",Rversion)
		if(grepl("<",Rversion)) matchtype <- "-lt"
		Rversion <- gsub("<","",Rversion)
		if(grepl("=",Rversion)) matchtype <- "-eq"
		Rversion <- gsub("=","",Rversion)
	
		num.version <- numeric(3)
		numout <- as.numeric(strsplit(Rversion, ".", fixed=TRUE)[[1]])
		num.version[1:length(numout)] <- numout
		num.version <- num.version * c(10^6,10^3,1)
		requiredvers <- num.version[1]+num.version[2]+num.version[3]
		
		# packages = either character string or list(runjags=">=1", coda=">1")
		if(class(packages)!="list" && class(packages)=="character"){
			pnames <- packages
			packages <- lapply(packages, function(x) return(">=0"))
			names(packages) <- pnames		
		}	
		packages <- lapply(packages, function(x){
			if(is.na(x) || length(x)==0 || !grepl("[[:digit:]]",x) || grepl("[[:alpha:]]",x) || grepl("<",x)){
				return(">=0")
			}else{
			 	return(x)
			}
			})
		packages <- packages[names(packages)!=""]
		
		art <- cat("#!/bin/bash
# For some reason which (and /usr/bin/which) seem to fail on private controllers at least, so use:
if [ ! -f '", Rpath, "' ]; then
	echo 0
else
	vers=`'", Rpath, "' --version | grep 'R version' | awk '{print $3}' | awk -F . '{print $1*1000000 + $2*1000 + $3 }'`
	if [ $vers ", matchtype, " ", format(requiredvers, scientific=FALSE), " ]; then", sep="", file=artfile, append=FALSE)
		
		# If we are looking for specific packages or need a (non default) art function to be assessed then:
		if(length(packages)>0 | !identical(body(eval(formals(setup.xgrid)$artfun)), body(artfun))){
				
			scriptfile <- tempfile()
			encodefile <- tempfile()
			dump("versionmatch", file=scriptfile)
			dump("packages", file=scriptfile, append=TRUE)
			dump("artfun", file=scriptfile, append=TRUE)
			
			cat("
installed <- installed.packages()[,'Version']
failed <- FALSE
for(i in 1:length(packages)){
	match <- which(names(installed) == names(packages)[i])
	if(length(match)==0 || !versionmatch(packages[[i]], installed[match])){
		cat('0\\n')
		failed <- TRUE
		break
	}
}
if(!failed){
	# Need to sink output when loading packages to make sure they can't produce any output to stdout/stderr:
	zz <- file(tempfile(), open='wt')
	sink(zz)
	sink(zz, type='message')
	success <- try({
	", if(length(packages)>0) paste("library(", paste(names(packages), collapse=", quietly=TRUE)\nlibrary("), ", quietly=TRUE)\n", sep=""), "
	})
	sink(type='message')
	sink()
	close(zz)
	unlink(zz)
	if(inherits(success, 'try-error')){
		writeLines('0')
	}else{
		artfun()
	}
}", sep="", file=scriptfile, append=TRUE)
		
			s <- system(paste('uuencode -m -o "', encodefile, '" "', scriptfile, '" script.R', sep=""), intern=TRUE, wait=TRUE)
	
			art <- cat("
		cat > rft.txt <<-EOA
	", sep="", file=artfile, append=TRUE)

		 	s <- system(paste('cat "', encodefile, '" >> "', artfile, '"', sep=""), intern=TRUE, wait=TRUE)

			cat("\nEOA
		uudecode -o script.R rft.txt
		R --slave -e \"source('script.R', echo=FALSE)\"
		rm rft.txt
		rm script.R", file=artfile, append=TRUE)
			# Clean up:
			unlink(encodefile)
			unlink(scriptfile)
						
		}else{
			cat("
		echo 1", file=artfile, append=TRUE)	
		}
		cat("
	else
		echo 0
	fi
fi
exit 0
", sep="", file=artfile, append=TRUE)
		
	}
	
	if(customrank){
		Sys.unsetenv("mgridcustomartname")
		system(paste('cp "', rankscript, '" "', artfile, '"', sep=""))
	}

	# Check the art file produced can be run:
	s <- system(paste('chmod 755 "', artfile, '"', sep=""), intern=TRUE)
	s <- as.numeric(system(paste('"', artfile, '"', sep=""), intern=TRUE))

	art <- paste(readLines(artfile), collapse="\n")
	unlink(artfile)
	
	if(s==0){
		# If the supplied artfun failes we have already had a warning so don't do it again:
		if(!artfails) warning(paste("The ART script created (with ", if(jagsrun) "JAGS" else "R", " version", if(!jagsrun) " and package", " dependencies) fails for this computer!", sep=""), call.=FALSE)
	}
	if(s!=1 & s!=0){
		stop(paste("An error occured while testing the ART script to enforce ", if(jagsrun) "JAGS" else "R", " version", if(!jagsrun) " and package", "dependencies", sep=""), call.=FALSE)
	}
	
	return(list(command=mgridcmd, method="xgrid", jobname=jobname, separate.tasks=separate.tasks, separate.jobs=separate.jobs, jagspath=jagspath, Rpath=Rpath, packages=packages, max.filesize=max.filesize, customart=art, submitandstop=submitandstop, mgridpath=mgridpath, hostname=hostname, password=password))
		
}


xgrid.retrieve <- function(jobinfo, wait, silent, cleanup, partialretrieve, jags=FALSE){
	
	if(class(jobinfo)=="list"){
		directory <- jobinfo$directory
		foldername <- jobinfo$jobname	
		jobnum <- jobinfo$jobid	
	}else{
		directory <- NA
		foldername <- jobinfo
		jobnum <- NA
	}
	if(class(foldername)!="character"){
		stop("The jobinfo must be supplied either as a character containing a path to the submitted folder, or the return value of an xgrid submit function", call.=FALSE)
	}
	
	# Find the relevant folder - first look for a matching folder in the working directory, and then go for the directory argument for an absolute path:
	if(!file.exists(foldername)){		
		if(!is.na(directory) && file.exists(directory)){
			foldername <- directory
		}else{
			stop("The folder specified by jobinfo does not exist", call.=FALSE)
		}
	}
	
	save.directory <- getwd()
	on.exit(setwd(save.directory))
	setwd(foldername)

	# JAGS stuff now no longer requires the .Rsave - this function is a horrible mess at the moment...
	if(jags){
		
		jobname <- jobinfo$jobname
		directory <- jobinfo$directory
		jobid <- jobinfo$jobid	
		
		n.sims <- 1
		separatejobs <- FALSE		
		
	}else{
		
		# Read the jobid and check it matches:
		if(!file.exists("info.Rsave")) stop("The required info.Rsave file was not found in the specified folder; ensure you have specified the path to a folder created by the xgrid submit functions", call.=FALSE)
	
		# To prevent warnings about no visible bindings:
		jagsrun = directory = jobname = jobid = threads = niters = iterations <- NA
		# In xgrid.submit and xgrid.jags.submit (minus threads,niters,iterations):
		# save(jagsrun, directory, jobname, jobid, threads, niters, iterations, file="info.Rsave")		
		load('info.Rsave')
	
		# Workingobj has an object in it called xg.jagsrun - ensure this matches with jags:
		if(jags != jagsrun) stop("Attempting to retrieve results from xgrid.submit using xgrid.results.jags (or vice versa) - ensure the results function type matches the submission function (JAGS or R)", call.=FALSE)
	
		# Note we ignore jobinfo$jobid - jobnum is saved in the workingobj.Rsave file:
		if(is.na(jobnum)){
			if(is.na(jobid)){
				stop("No valid job ID was supplied or found in the xgrid job directory", call.=FALSE)
			}else{
				jobnum <- jobid
			}
		}
	
		if(length(jobnum) > 1) separatejobs <- TRUE else separatejobs <- FALSE
		n.sims <- length(jobnum)
		
	}
	
		
	if(separatejobs){
		
		status <- character(n.sims)
		done <- replicate(n.sims, FALSE)
		for(s in 1:n.sims){
			statusout <- system(paste('xgrid -job attributes -id ', jobnum[s], sep=''), intern=TRUE)
			if(paste(statusout, collapse='')=="{    error = InvalidJobIdentifier;}") stop("One of the jobs specified is not on xgrid.  This can sometimes occur when a job is being initialised - you could try again in a minute or so...", call.=FALSE)
			tstatus <- statusout[grep('jobStatus', statusout)]
			status[s] <- gsub('[[:space:]]', '', gsub(';', '', gsub('jobStatus = ', '', tstatus)))
			if(any(status[s]==c('Finished', 'Canceled', 'Failed'))) done[s] <- TRUE
		}
		
		if(!wait & !all(done) & !partialretrieve){
			if(!silent){
				swcat('Jobs not finished.  Statuses are "', paste(status, collapse=', '), '"\n\nThe job outputs (if any) are shown below:\n', sep='')
				for(s in 1:n.sims){
					system(paste('xgrid -job results -id ', jobnum[s], sep=''))
					swcat('\n')
				}
			}
			stop(paste('Jobs not finished.  Statuses are "', paste(status, collapse=', '), '"', sep=''))
		}
		
		if(wait & !all(done)){			
			swcat('Job statuses at ', format(Sys.time(), "%a %b %d %H:%M:%S"), ' were "', paste(status, collapse=', '), '".  Waiting for jobs to complete....\n', sep='')
			
			for(s in 1:n.sims){
				res <- waitforxgridjob(j=jobnum[s])
			}
			for(s in 1:n.sims){
				statusout <- system(paste('xgrid -job attributes -id ', jobnum[s], sep=''), intern=TRUE)
				if(paste(statusout, collapse='')=="{    error = InvalidJobIdentifier;}") stop("One or more of the jobs has been deleted from xgrid", call.=FALSE)
				tstatus <- statusout[grep('jobStatus', statusout)]
				status[s] <- gsub('[[:space:]]', '', gsub(';', '', gsub('jobStatus = ', '', tstatus)))
				if(any(status[s]==c('Finished', 'Canceled', 'Failed'))) done[s] <- TRUE				
			}
		}

		if(all(status=='Finished')){
			swcat('The xgrid jobs have finished\n')
		}else{
			swcat('The xgrid jobs are showing the statuses "', paste(status, collapse=', '), '"\n', sep="")
			silent <- FALSE
		}
		
		xgridoutput <- vector('list', length=n.sims)
		for(s in 1:n.sims){
			xgridoutput[[s]] <- c(paste('\n', if(jags) 'Chain' else 'Task', ' ', s, ':', sep=''), system(paste('xgrid -job results -id ', jobnum[s], ' -out "', directory, if(jags) '/sim.', if(jags) s, '"', sep=''), intern=TRUE))
			if(length(xgridoutput[[s]])==0 & jags) stop(paste("The job produced no output for chain ", s, "; ensure that the jagspath supplied is accurate", sep=''), call.=FALSE)
		}

		if(!silent)	cat('Job was successfully retreived from xgrid\n\nThe xgrid output is displayed below:\n', unlist(xgridoutput), sep='\n')
		if(cleanup){
			for(s in 1:n.sims){
				xgriddeleteout <- system(paste('xgrid -job delete -id ', jobnum[s], sep=''), intern=TRUE)
				if(paste(xgriddeleteout, collapse='')!='{}') warning(paste('Possible error deleting xgrid job number ', jobnum[s], ' - please check this manually', sep=''))
			}
		}
		
	}else{
		
		statusout <- system(paste('xgrid -job attributes -id ', jobnum, sep=''), intern=TRUE)
		if(paste(statusout, collapse='')=="{    error = InvalidJobIdentifier;}") stop("The job specified is not on xgrid.  This can sometimes occur when a job is being initialised - you could try again in a minute or so...", call.=FALSE)
		status <- statusout[grep('jobStatus', statusout)]
		status <- gsub('[[:space:]]', '', gsub(';', '', gsub('jobStatus = ', '', status)))
		
		if(!wait & !status=='Finished' & !partialretrieve){
			if(!silent){
				swcat('Job not finished.  Status is "', status, '"\n\nThe job output (if any) is shown below:\n')
				system(paste('xgrid -job results -id ', jobnum, sep=''))
			}
			stop(paste('Job not finished.  Status is "', status, '"', sep=''), call.=FALSE)
			done <- FALSE
		}else{
			done <- TRUE
		}
		
		if(!any(status==c('Finished', 'Canceled', 'Failed')) & wait){
			swcat('Job status at ', format(Sys.time(), "%a %b %d %H:%M:%S"), ' was "', status, '".  Waiting for job to complete....\n', sep='')
		
			res <- waitforxgridjob(jobnum)
									
			statusout <- system(paste('xgrid -job attributes -id ', jobnum, sep=''), intern=TRUE)
			if(paste(statusout, collapse='')=="{    error = InvalidJobIdentifier;}") stop("The job has been deleted from xgrid")
			status <- statusout[grep('jobStatus', statusout)]
			status <- gsub('[[:space:]]', '', gsub(';', '', gsub('jobStatus = ', '', status)))
			if(any(status==c('Finished', 'Canceled', 'Failed'))) done <- TRUE
						
		}
	
		if(status=='Finished'){
			swcat('The xgrid job has finished\n')
		}else{
			swcat('The xgrid job is showing the status "', status, '"\n', sep="")
			silent <- FALSE
		}
		xgridoutput <- system(paste('xgrid -job results -id ', jobnum, ' -out "', directory, '"', sep=''), intern=TRUE)
		
		swcat("Job was successfully retreived from xgrid\n")
		if(!silent){
			# Note not using swcat for xgrid output!
			if(length(xgridoutput)==0) swcat('\nThe xgrid job has not produced any output to screen\n') else cat('\nThe xgrid output is displayed below:\n', xgridoutput, sep='\n')
		}
		
		if(cleanup){
			xgriddeleteout <- system(paste('xgrid -job delete -id ', jobnum, sep=''), intern=TRUE)
			if(paste(xgriddeleteout, collapse='')!='{}') warning(paste('Possible error deleting xgrid job number ', jobnum, ' - please check this manually', sep=''))
		}
		
	}
	
	return(list(done=done, jobname=foldername, directory=directory, xgridoutput=xgridoutput))
}

waitforxgridjob <- function(j){
	done <- FALSE
	tryCatch({
		flush.console()
		if(.Platform$GUI=="AQUA"){
			s <- system(paste("xgrid -job wait -id ", j, " 2>&1 > waitfile.txt", sep=""), intern=FALSE, wait=FALSE)
			res <- tailf("waitfile.txt", stop.text="jobStatus = ", print=FALSE, return=TRUE, refresh=2)$text		
		}else{
			res <- system(paste("xgrid -job wait -id ", j, sep=""), intern=TRUE, wait=TRUE)	
		}
		# Sending interrupt to system calls doesn't trigger a stop, so check the xgrid call output something when it finished:
		if(grepl("jobStatus = ", paste(res,collapse="\n"))) done <- TRUE

	}, finally={

		if(!done) stop("User cancelled waiting for xgrid job to finish - the job has not been deleted", call.=FALSE)
		
		return(res)	
	})
	
}

