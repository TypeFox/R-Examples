ask <- function (prompt="?", type="logical", bounds=c(-Inf, Inf), na.allow=FALSE){
	prompt <- paste(prompt, '  ', sep='')
	
	rettypes <- c("logical", "numeric", "character", "integer")
	type <- pmatch(type,rettypes)	
	if(is.na(type)){
		stop(paste("Unrecognised return type \'", type, "\'", sep=""))
	}
	type <- rettypes[type]
	
	if(type=="logical"){
		repeat{
			result <- readline(prompt = prompt)
			if(na.allow==TRUE && (result=="NA" | result=="na" | result=="Na")){
				return(NA)
			}
			if((result == "T") | (result == "TRUE") | (result == "True") | (result == "true") | (result == "t") | (result == "y") | (result == "YES") | (result == "yes") | (result == "Yes") | (result == "Y")){
				return(TRUE)
			}
			if((result == "F") | (result == "FALSE") | (result == "False") | (result == "false") | (result == "f") | (result == "n") | (result == "NO") | (result == "no") | (result == "No") | (result == "N")){
				return(FALSE)
			}
			swcat("ERROR:  Please enter 'y' or 'n'\n")
		}
	}
	
	if(type=="numeric"){
		repeat{
			result <- readline(prompt = prompt)
			if(na.allow==TRUE && (result=="NA" | result=="na" | result=="Na")){
				return(NA)
			}
			suppressWarnings(result <- as.numeric(result))
			if(is.na(result)){
				swcat("ERROR:  Please enter a number\n")
			}else{
				if((result > bounds[2]) | (result < bounds[1])){
					swcat("ERROR:  Please enter a number between ", bounds[1], " and ", bounds[2], "\n", sep="")
				}else{
					return(result)
				}
			}
		}
	}
	
	if(type=="character"){
		return(as.character(readline(prompt=prompt)))
	}
	
	if(type=="integer"){
		repeat{
			result <- readline(prompt=prompt)
			if(na.allow==TRUE && (result=="NA" | result=="na" | result=="Na")){
				return(NA)
			}
			suppressWarnings(result <- as.numeric(result))
			if(is.na(result) | (as.integer(result) != result)){
				swcat("ERROR:  Please enter a whole number\n")
			}else{
				if((result > bounds[2]) | (result < bounds[1])){
					swcat("ERROR:  Please enter a whole number between ", bounds[1], " and ", bounds[2], "\n", sep="")
				}else{
					return(result)
				}
			}
		}
	}			
}

findjags <- function(ostype = .Platform$OS.type, look_in = NA, ...){
	
	if(is.na(look_in)) look_in <- if(ostype=="windows") c("/Program Files/","/Windows/Program Files/","C:/Program Files/","C:/Windows/Program Files/","/") else NULL
		
	if(! ostype %in% c("unix","windows")) stop(paste("Unrecognised OS type '", ostype, "'.  Use either 'unix' or 'windows'", sep=""))
	
	# This is a bit of a mess now the option is deprecated...
	if(!exists('from.variable')) from.variable <- ""
	if(is.na(from.variable) | is.null(from.variable) | from.variable=="" | length(from.variable)==0){
		from.variable <- character(0)
	}else{
		if(length(find(from.variable))==0){
			from.variable <- character(0)
		}else{
			from.variable <- get(from.variable)
		}		
	}
	if(length(from.variable)!=0){
		warning("The '.jagspath' option is deprecated - use the runjags.options() to set the default path to JAGS")
	}

	if(ostype=="unix"){
		
		# Test standard 'jags' as well as variable passed and some common install locations:
		suppressWarnings({
			paths <- c(from.variable, "jags", "/opt/local/bin/jags", "/opt/local/sbin/jags", "/usr/texbin/jags", "/usr/bin/jags", "/bin/jags", "/usr/sbin/jags", "/sbin/jags", "/usr/local/bin/jags", "/usr/X11/bin/jags")
			for(i in 1:length(paths)){
				jagspath <- system(paste('which ', paths[i], ' 2>&1', sep=""), intern=TRUE)
				if(length(jagspath)!=0) return(jagspath)
			}
		})
		
	}
	
	if(ostype=="windows"){
		
		# First check if the environmental variable is set (by rjags), then this will work:
		suppressWarnings(s <- try(system('where jags',intern=TRUE),silent=TRUE))
		if(is.null(attr(s,'status'))){
			s <- gsub('jags.bat','jags-terminal.exe',s)
			s <- gsub('\\','/',s,fixed=TRUE)
			return(s)
		}
		
		if(class(look_in)!="character") stop("The look_in argument supplied to findjags must be a character vector of paths")
		look_in <- paste(look_in,"/",sep="")
		look_in <- gsub("//","/",look_in,fixed=TRUE)
		suppressWarnings(paths <- paste(unlist(lapply(look_in, function(x) return(paste(x,list.files(x),sep="")))),"/",sep=""))
		possible <- unique(paths[grepl("JAGS",paths)|grepl("jags",paths)])
		
		if(Sys.getenv("JAGS_HOME")!="") posspaths <- Sys.getenv("JAGS_HOME") else posspaths <- character(0)

		if(length(possible)>0) for(i in 1:length(possible)){
		binpaths <- list.files(possible[i],recursive=TRUE)
		posspaths <- c(posspaths, paste(possible[i],binpaths[grepl("jags-terminal.exe",binpaths)],sep=""))
		}

		posspaths <- c(unique(posspaths), from.variable)
		versions <- sapply(posspaths, function(x){
			pathsegs <- strsplit(x,"/")[[1]]
			withnum <- grepl("[0-9]",pathsegs)
			return(getvers(pathsegs[withnum][1]))
		})
		# If versions doesn't contain a number, getvers will return 0

		x64s <- grepl('x86',posspaths)

		if(any(x64s)){
			versions <- versions[x64s]
			posspaths <- posspaths[x64s]
		}
		
		if(length(na.omit(versions))>0){
			return(posspaths[which(versions==max(versions))[1]]) 
		}
		
		##  jags.bat spawns lots of horrible windows.  jags-terminal.exe requires environmental variables - set by testjags()
	}
	
	# JAGS executable wasn't found - see if rjags is installed:
	
	# If rjags is loaded no warning:
	if(! any(.packages(FALSE)=='rjags')){
		# If rjags is installed but not loaded suggest loading it:
		
		if(any(.packages(TRUE)=='rjags')){
			warning("JAGS was not found in the common install locations on your system; try loading the rjags package first or alternatively provide a path to the executable as the jags argument or using runjags.options()")
		}else{
			warning("JAGS was not found in the common install locations on your system; please provide a path to the executable as the jags argument or using runjags.options()")
		}
	}
	
	return("JAGS not found")
	
}


new_unique <- function(name=NA, suffix="", ask=FALSE, prompt="A file or directory with this name already exists.  Overwrite?", touch=FALSE, type='file'){
	
	if(suffix!=""){
		if(strsplit(suffix, split="")[[1]][1]!=".") suffix <- paste(".", suffix, sep="")
	}
	
	if(!any(type==c('file', 'folder', 'directory', 'dir', 'f', 'd'))) stop('The type specified must be either file or directory')
	
	if(type=='f') type <- 'file'
	if(any(type==c('folder', 'dir', 'd'))) type <- 'directory'
	
	for(i in 1:length(name)){
		if(is.na(name[i])){
			alphanumeric <- c(as.numeric(0:9), "a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z", toupper(c("a", "b", "c", "d", "e", "f", "g", "h", "i", "j", "k", "l", "m", "n", "o", "p", "q", "r", "s", "t", "u", "v", "w", "x", "y", "z")))
			name[i] <- paste(sample(alphanumeric, 9, replace=TRUE), sep="", collapse="")
		}
	}
	
	temp <- paste(name, sep="", collapse="")
	exists <- file.exists(paste(temp, suffix, sep=""))
	if(exists==TRUE){
		path.ok <- FALSE
		counter <- 1
		if(ask==TRUE){
			path.ok <- ask(paste("\'", temp, suffix, "\'.  ", prompt, "  ", sep=""), type="logical")
			if(path.ok==TRUE){
				unlink(paste(temp, suffix, sep=""), recursive = TRUE)
			}
		}
		while(path.ok == FALSE){
			temptemp <- paste(temp, "_", counter, "", sep="")
			exists <- file.exists(paste(temptemp, suffix, sep=""))
			if(exists==TRUE){
				counter <- counter + 1
			}else{
				path.ok <- TRUE
				temp <- temptemp
				break
			}
		}
	}
	
	if(type=='file'){
		suppressWarnings(try(file.create(paste(temp, suffix, sep="")), silent=TRUE))
	}else{
		suppressWarnings(try(dir.create(paste(temp, suffix, sep="")), silent=TRUE))
	}
	permissions <- file.exists(paste(temp, suffix, sep=""))
	if(permissions==FALSE){
		swcat("Error:  Directory not writable\n")
		return("Directory not writable")
	}else{
		if(!touch) unlink(paste(temp, suffix, sep=""), recursive = TRUE)
	}
	if(!touch) suppressWarnings(try(backupforspaces <- file.remove(paste(temp, suffix, sep="")), silent=TRUE))
	return(paste(temp, suffix, sep=""))
}


timestring <- function(time1, time2=NA, units=NA, show.units=TRUE){
	
	time <- na.omit(c(time1, time2))
	
	if(length(time)==2){
		time <- as.integer(difftime(time[2], time[1], units="secs")*10)
		time <- time/10
	}else{	
		if(length(time)==1){
			time <- as.integer(time*10, units="secs")/10
		}else{
			stop("Input variables incorrectly specified")
		}
	}
	
	secs <- time
	mins <- round(time / (60), digits=1)
	hours <- round(time / (60*60), digits=1)
	days <- round(time / (60*60*24), digits=1)
	weeks <- round(time / (60*60*24*7), digits=1)	
	years <- round(time / (60*60*24*7*52), digits=1)
	
	if(!is.na(units)){
		if(units=="s" | units=="secs" | units=="seconds"){
			if(show.units==FALSE){
				return(secs)
			}else{
				return(paste(secs, " seconds", sep=""))
			}
		}
		if(units=="m" | units=="minutes" | units=="mins"){
			if(show.units==FALSE){
				return(mins)
			}else{
				return(paste(mins, " minutes", sep=""))
			}
		}
		if(units=="h" | units=="hours"){
			if(show.units==FALSE){
				return(hours)
			}else{
				return(paste(hours, " hours", sep=""))
			}
		}
		if(units=="d" | units=="days"){
			if(show.units==FALSE){
				return(days)
			}else{
				return(paste(days, " days", sep=""))
			}
		}
		if(units=="w" | units=="weeks"){
			if(show.units==FALSE){
				return(weeks)
			}else{
				return(paste(weeks, " weeks", sep=""))
			}
		}
		if(units=="y" | units=="years"){
			if(show.units==FALSE){
				return(years)
			}else{
				return(paste(years, " years", sep=""))
			}
		}
		swcat("Error:  Unrecognised unit type '", units, "'\n", sep="")
	}
		
	if(secs < 60){
		if(show.units==FALSE){
			return(secs)
		}else{
			return(paste(secs, " seconds", sep=""))
		}
	}
	if(mins < 60){
		if(show.units==FALSE){
			return(mins)
		}else{
			return(paste(mins, " minutes", sep=""))
		}
	}
	if(hours < 24){
		if(show.units==FALSE){
			return(hours)
		}else{
			return(paste(hours, " hours", sep=""))
		}
	}
	if(days < 7){
		if(show.units==FALSE){
			return(days)
		}else{
			return(paste(days, " days", sep=""))
		}
	}
	if(weeks < 52){
		if(show.units==FALSE){
			return(weeks)
		}else{
			return(paste(weeks, " weeks", sep=""))
		}
	}
	
	if(show.units==FALSE){
		return(years)
	}else{
		return(paste(years, " years", sep=""))
	}
	
}

testjags <- function(jags=runjags.getOption('jagspath'), silent=FALSE){
	
	# Evaluate jags outside a suppresswarnings:
	jags <- jags
	
	tempfile <- new_unique(name=tempfile("test"), suffix=".cmd", touch=TRUE, type='file')
	write("exit", file=tempfile)
	s.info <- Sys.info()
	p.info <- .Platform
	
	os <- p.info$OS.type
	username <- as.character(s.info["user"])
	rversion <- R.version$version
	gui <- p.info$GUI
	p.type <- p.info$pkgType
	
	libpaths <- NULL
	if(os=="windows"){
		
		if(file.exists(paste(jags, 'bin', .Platform$file.sep, 'jags-terminal.exe', sep=''))) jags <- paste(jags, 'bin/jags-terminal.exe', sep='')
		if(file.exists(paste(jags, .Platform$file.sep, 'bin', .Platform$file.sep, 'jags-terminal.exe', sep=''))) jags <- paste(jags, '/bin/jags-terminal.exe', sep='')
		if(file.exists(paste(jags, 'jags-terminal.exe', sep=''))) jags <- paste(jags, 'jags-terminal.exe', sep='')
		if(file.exists(paste(jags, .Platform$file.sep, 'jags-terminal.exe', sep=''))) jags <- paste(jags, '/jags-terminal.exe', sep='')
				
		path <- strsplit(jags, split=paste(.Platform$file.sep, 'bin', sep=''), fixed=TRUE)

		firstpath <- path[[1]][1:length(path[[1]])-1]

		binpath <- paste(firstpath, .Platform$file.sep, 'bin;', sep='')
		libpath <- paste(firstpath, .Platform$file.sep, 'modules;', sep='')
		
		libpaths <- list(PATH=binpath, LTDL_LIBRARY_PATH=libpath)
	}
	
	suppressWarnings(returnval <- try(system2(jags, args=tempfile, stdout=TRUE, stderr=TRUE), silent=TRUE))
	unlink(tempfile)
	
	if(class(returnval)=="try-error"){		
		success <- 0
	}else{
		if(is.null(attributes(returnval))){
			success <- 1			
		}else{
			success <- -1
		}
	}
	
	if(length(returnval)==0){
		success <- 0
#		warning(paste("Running the command ", jags, " returned no output - ensure the path is correct and possibly try reinstalling jags\n", sep=""))
	}
	if(!any(grepl("Welcome",returnval))){
		success <- 0
#		warning(paste("Running the command ", jags, " did not return the expected welcome message - ensure the path is correct and possibly try reinstalling jags\n", sep=""))
	}
	
	if(success){
    	rightstring <- which(grepl("Welcome",returnval))[1]
    	if(is.na(rightstring)){  # Will be NA if which is length 0 as selected first element above
    		version <- 'unknown'
    		num.version <- NA
    	}else{
    		version <- strsplit(returnval[rightstring], split=" ", fixed=TRUE)[[1]][4]
			versioncom <- version
			if(grepl('(', versioncom, fixed=TRUE)){
				versioncom <- strsplit(versioncom, split='(', fixed=TRUE)[[1]][1]
			}
			versioncom <- sub('.',';',versioncom,fixed=TRUE)
			versioncom <- gsub('.','',versioncom,fixed=TRUE)
			num.version <- as.numeric(gsub(';','.',versioncom,fixed=TRUE))
			if(is.na(version)) version <- 'unknown'
    	}
    }
    
	if(any(.packages(TRUE)=="rjags")){
		rjags.avail <- TRUE
		rjags.major <- packageVersion('rjags')$major
		rjags.version <- packageDescription('rjags', fields='Version')
	}else{
		rjags.avail <- FALSE
		rjags.major <- 0
		rjags.version <- ""
	}
	
	if(!silent)
		swcat("You are using ", rversion, " on a ", os, " machine, with the ", gui, " GUI\n", sep="")

	
	if(jags=="JAGS not found") jags <- "findjags()"
		
	if(success){
		if(success==-1){
			if(!silent) suppressWarnings(swcat("JAGS version ", version, " found successfully using the command '", jags, "', but returned the status code '", attributes(returnval)$status, "' - this may indicate a compatibility issue, procede with caution\n", sep=""))
		}else{
			if(!silent) suppressWarnings(swcat("JAGS version ", version, " found successfully using the command '", jags, "'\n", sep=""))
		}
		if(!is.na(num.version) & num.version<1){
			warning(paste("The version of JAGS currently installed on your system is no longer supported.  Please update JAGS from http://mcmc-jags.sourceforge.net\n"))
			jags.avail <- FALSE
		}else{
			jags.avail <- TRUE
		}
		jagsfound <- TRUE
	}else{
		if(rjags.avail){
			if(!silent) suppressWarnings(swcat("JAGS was not found on your system using the command '", jags, "'.  Please ensure that the command is correct and that the latest version of JAGS is installed from http://mcmc-jags.sourceforge.net\n", sep=""))
			jags.avail <- TRUE
			# If it's just rjags assume the version number is high enough:
			num.version <- Inf
		}else{
			if(!silent) suppressWarnings(swcat("JAGS was not found on your system using the command '", jags, "'.  Please ensure that the command is correct.\n", sep=""))
			jags.avail <- FALSE
			num.version <- "none found"			
		}
		if(gui %in% c('RStudio', 'AQUA') && !silent)
			swcat("Note that as of OS X Yosemite, the FULL path must be provided to JAGS as the $PATH global variable is not passed to processes started from within a GUI application\n", sep="")
		
		jagsfound <- FALSE
	}
	
	if(rjags.avail && num.version < Inf){  # If rjags and JAGS are both findable independently
		if(success!=0 && rjags.major!=floor(num.version)){   # If JAGS is installed and rjags major != JAGS major
			if(!silent)
				swcat("The rjags package version ", rjags.major, ".x is installed, but this is not compatible with your installed version of JAGS (version ", floor(num.version), ".x)\n",sep="")

			if(rjags.major > num.version){
				msg <- paste("You should update to JAGS version ", rjags.major, ".x from https://sourceforge.net/projects/mcmc-jags/files/JAGS/", sep="")
				
			}else{
				msg <- paste("You should update the rjags package to version ", floor(num.version), ".x - if the version on CRAN is not currently up to date, try downloading from https://sourceforge.net/projects/mcmc-jags/files/rjags/ instead", sep="")
			}
			if(!silent)
				swcat(msg,"\n")

			if(!runjagsprivate$warned_version_mismatch){
				warning(paste(msg, "\n(This warning is given once per R session)"), call. =FALSE)
				runjagsprivate$warned_version_mismatch <- TRUE
			}
			
		}else if(!silent){
			swcat("The rjags package is installed\n",sep="")
		}
	}else if(rjags.avail && num.version == Inf && !silent){
		swcat("The rjags package version ", rjags.major, ".x is installed (compatibility with the major version of JAGS cannot be verified)\n",sep="")
	}else if(!rjags.avail && !silent){
		swcat("The rjags package is not installed\n",sep="")
	}
		
	invisible(list("os"=os, "JAGS.available"=jags.avail, "JAGS.found"=jagsfound, "rjags.found"=rjags.avail, "rjags.version"=rjags.version, "rjags.major"=rjags.major, "JAGS.path"=jags, "JAGS.version"=version, "JAGS.major"=floor(num.version), "R.version"=rversion, "R.GUI"=gui, "R.package.type"=p.type, "username"=username, libpaths=libpaths))
}

findJAGS <- findjags
testJAGS <- testjags



runjagsprivate <- new.env()
# Use 'expression' for functions to avoid having to evaluate before the package is fully loaded:
assign("defaultoptions",list(jagspath=expression(findjags()),method=expression(if('rjags' %in% .packages(TRUE)){'rjags'}else{if(Sys.info()['user']=='nobody') 'simple' else 'interruptible'}), tempdir=TRUE, plot.layout=c(2,2), new.windows=TRUE, modules="", factories="", bg.alert='beep', linenumbers=TRUE, inits.warning=TRUE, rng.warning=TRUE, summary.warning=TRUE, blockcombine.warning=TRUE, blockignore.warning=TRUE, tempdir.warning=FALSE, nodata.warning=TRUE, adapt.incomplete='warning', repeatable.methods=FALSE, silent.jags=FALSE, silent.runjags=FALSE, predraw.plots=FALSE, force.summary=FALSE, mode.continuous=expression('modeest' %in% .packages(TRUE)), timeout.import=30, partial.import=FALSE, keep.crashed.files=TRUE, full.cleanup=FALSE, debug=FALSE), envir=runjagsprivate)
assign("options",runjagsprivate$defaultoptions,envir=runjagsprivate)
assign("rjagsmethod",c('rjags','rjparallel'),envir=runjagsprivate)
assign("bgmethod",c('background','bgparallel'),envir=runjagsprivate)
assign("parallelmethod",c('parallel','bgparallel','snow','rjparallel','xgrid'),envir=runjagsprivate)
assign("runjagsversion", "notset", envir=runjagsprivate)
assign("simfolders", character(0), envir=runjagsprivate)
assign("failedsimfolders", character(0), envir=runjagsprivate)
assign("defaultsummarypars", list(vars=NA, mutate=NULL, psrf.target = 1.05, normalise.mcmc = TRUE, modeest.opts=list(), confidence=c(0.95), autocorr.lags=c(10), custom=NULL, silent.jags=expression(runjags.getOption('silent.jags')), plots=FALSE, plot.type=c('trace','ecdf','histogram','autocorr','key','crosscorr'), col=NA, summary.iters=10000, trace.iters=1000, separate.chains=FALSE, trace.options=list(), density.options=list(), histogram.options=list(), ecdfplot.options=list(), acplot.options=list()), envir=runjagsprivate)
assign("minjagsmajor", 3, envir=runjagsprivate)
assign("maxjagsmajor", 4, envir=runjagsprivate)
assign("warned_version_mismatch", FALSE, envir=runjagsprivate)

	# runjags.getOption is not available at compile time so has to be expression, but it's OK as it is eval()ed when getting defaultsumpars
getdefaultsummarypars <- function(){
	x <- runjagsprivate$defaultsummarypars
	x$silent.jags <- eval(x$silent.jags)
	return(x)
}
assign("defaultmethodoptions", list(n.sims=NA, cl=NA, remote.jags="*//*usefindjags*//*", rjags=NA, by=NA, progress.bar=expression(unlist(options('jags.pb'))), jags=expression(runjags.getOption('jagspath')), silent.jags=expression(runjags.getOption('silent.jags')), jags.refresh=0.1, batch.jags=expression(runjags.getOption('silent.jags'))), envir=runjagsprivate)
	# jags.pb is set by runjags - but its ok as it is eval()uated before use
	# is.na(cl/by) is what is expected if they are to be set at run time
getdefaultmethodoptions <- function(){
	x <- runjagsprivate$defaultmethodoptions
	x$progress.bar <- eval(x$progress.bar)
	x$silent.jags <- eval(x$silent.jags)
	x$jags <- eval(x$jags)
	x$batch.jags <- eval(x$batch.jags)
	if(is.null(x$progress.bar))   # Will be NULL if rjags isn't installed
    x$progress.bar <- 'text'
	return(x)
}

runjags.options <- function(...){
	opts <- list(...)

	# For backwards compatibility:
	if(any(names(opts)=='newwindows')){
		opts$new.windows <- opts$newwindows
		opts$newwindows <- NULL
	}

	if(length(opts)>0){
		options <- runjagsprivate$options
		recognised <- pmatch(names(opts), names(options))
		if(any(is.na(recognised))){
			warning(paste("Igoring unmatched or ambiguous option(s): ", paste(names(opts)[is.na(recognised)],collapse=", ")))
      opts <- opts[!is.na(recognised)]
		}
		optnames <- names(options)[recognised[!is.na(recognised)]]
		if(length(optnames)>0) for(i in 1:length(optnames)){
			options[optnames[i]] <- opts[[i]]
		}
		assign("options",options,envir=runjagsprivate)
	}
	
	# Some checks for valid option settings:
	if((!is.numeric(runjags.getOption('plot.layout')) || length(runjags.getOption('plot.layout'))!=2))
		stop('The "plot.layout" option must be a numeric vector of length 2')
	
	runjagsprivate$options$adapt.incomplete <- c('error','warning','silent')[pmatch(runjagsprivate$options$adapt.incomplete,c('error','warning','silent'))]
	if(is.na(runjagsprivate$options$adapt.incomplete))
		stop('The "adapt.incomplete" option must be either error, warning or silent')
	
	invisible(runjagsprivate$options)
}

runjags.getOption <- function(name){
	if(length(name)!=1) stop("Only 1 option can be retrieved at a time")
	opt <- pmatch(name,names(runjagsprivate$options))
	if(is.na(opt)) stop(paste("Unmatched or ambiguous option '", name, "'", sep=""))
	# Use eval as some defaults are put in using 'expression' to avoid evaluating at load time:
	return(eval(runjagsprivate$options[[opt]]))
}

runJAGS.getOption <- runjags.getOption


failedjags <- new.env()

assign("model", NA, envir=failedjags)
assign("data", NA, envir=failedjags)
assign("inits", NA, envir=failedjags)
assign("output", NA, envir=failedjags)
assign("end.state", NA, envir=failedjags)

