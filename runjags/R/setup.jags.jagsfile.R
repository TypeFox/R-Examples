setup.jags <- function(model, monitor = stop("No monitored variables supplied"), data=NA,  n.chains=2, inits = replicate(n.chains, NA), modules=c(""), factories=c(""), response=NA, fitted =NA, residual=NA, jags = runjags.getOption('jagspath'), method="simple", mutate=NA){
	
	# Reset failedjags stuff:
	failedjags$model <- NA
	failedjags$data <- NA
	failedjags$inits <- NA
	failedjags$output <- NA
	failedjags$end.state <- NA
	
	# We may be passed some unevaluated function arguments so evaluate everything here:
	argnames <- names(formals(setup.jags))
	for(i in 1:length(argnames)){
		success <- try(assign(argnames[i], eval(get(argnames[i]))), silent=TRUE)		
		if(inherits(success, 'try-error')){
			stop(paste("object '", strsplit(as.character(success),split="'",fixed=TRUE)[[1]][2], "' not found", sep=""), call.=FALSE)
		}
	}
	
  # Convert specified character to lists - or '' if none:
	factories <- checkmodfact(factories, 'factory')
	modules <- checkmodfact(modules, 'module')
	
	if(inherits(model, "runjagsmodel")) class(model) <- "character"
	if(inherits(data, "runjagsdata")) class(data) <- "character"
	if(inherits(inits, "runjagsinits")) class(inits) <- "character"
	
	if(identical(inits, list()) || is.null(inits))
		inits <- NA
	
	jags.status <- testjags(jags, silent=TRUE)
	if(jags.status$JAGS.available==FALSE){
		if(jags.status$os=="windows"){
			# Try it again - sometimes this seems to clear it up:
			Sys.sleep(0.2)
			jags.status <- testjags(jags, silent=TRUE)
		}		
		jags <- jags.status$JAGS.path
		
		if(jags.status$JAGS.available==FALSE){			
			swcat("Unable to call JAGS using '", jags, "' - try specifying the path to the JAGS binary as the jags argument, or installing the rjags package.  Use the testjags() function for more detailed diagnostics.\n", sep="")
			stop("Unable to call JAGS", call.=FALSE)
		}
	}
	jags <- jags.status$JAGS.path
		
	if(!jags.status$JAGS.found && ! method%in%c("snow",runjagsprivate$rjagsmethod)){
		swcat("Unable to call JAGS using '", jags, "' - try specifying the path to the JAGS binary as the jags argument, or using the rjags method.  Use the testjags() function for more detailed diagnostics.\n", sep="")
		stop("Unable to call JAGS", call.=FALSE)
	}

	if(method %in% runjagsprivate$rjagsmethod)
		rjagsmethod <- TRUE else rjagsmethod <- FALSE
	if(rjagsmethod)
		loadandcheckrjags()
	
	if(class(monitor)!="character" | all(is.na(monitor))){
		stop("Monitored variable(s) must be provided in the form of a character vector")
	}
	monitor <- checkvalidmonitorname(monitor)		

	if(any(tolower(monitor)=="deviance") && (any(tolower(monitor)=="full.pd") || any(tolower(monitor)=="pd"))){
		monitor <- c(monitor, "dic")
	}
	if(any(tolower(monitor)=="deviance") && any(tolower(monitor)=="popt")){
		monitor <- c(monitor, "ped")
	}
	
	monitor[monitor=="DIC"] <- "dic"
	monitor[monitor=="PED"] <- "ped"
	monitor[monitor=="pD"] <- "pd"
	monitor[monitor=="full.pD"] <- "full.pd"
	monitor[monitor=="pOpt"] <- "popt"	

	monitor[monitor=="pD.i"] <- "pd.i"
	monitor <- na.omit(monitor[monitor!=""])
	monitor <- (unique(monitor))
	
	if(any(monitor=='pd.i')){
		warning('Ignoring deprecated pd.i monitor', call.=FALSE)
		monitor <- monitor[monitor!='pd.i']
	}
	
	if(any(c("popt", "pd", "full.pd", "deviance", "dic", "ped") %in% monitor)){
		if(identical(modules, ''))
			modules <- list(c("dic","TRUE"))
		if(!any(sapply(modules,function(x) return(x[1]))=='dic'))
			modules <- c(modules, list(c("dic","TRUE")))
	}
	monitor <- (unique(monitor))
	
	if(class(model)!="character" | length(model)!=1){
		stop("The model must be provided in the form of a character string")
	}
	
	# Find references to functions in the runjags module (unless we have already been told to use the standalone runjags module):
	if(!any(modules=="paretoprior")){
		# Find any matching functions used (not in comments, and definitely used as a function not a variable):
		fs <- c("par1","par2","par3","par4","lomax","mouch","genpar","halfcauchy")

		fs <- apply(expand.grid(c("~","<-"),c("d","p","q","logdensity."),fs,"("),1,paste,collapse="")
		# Get rid of commented lines and remove all spaces:
		nohashstring <- paste(lapply(strsplit(model, "[\n\r]")[[1]], function(x) gsub("#.*", "", x)), collapse="\n")
		nohashstring <- gsub('[[:space:]]','',nohashstring)

		# Find any matches and add runjags to the modules if rjags method:
		if(any(sapply(fs,function(x) return(grepl(x,nohashstring,fixed=TRUE))))){
			if(rjagsmethod){
				if(identical(modules, ''))
					modules <- list(c("runjags","TRUE"))				
				if(!any(sapply(modules,function(x) return(x[1]))%in%c('runjags','paretoprior')))
					modules <- c(modules, list(c("runjags","TRUE")))
			}else{
				if(!any(sapply(modules,function(x) return(x[1]))%in%c('paretoprior')))
					warning("Pareto family functions provided by the runjags module are only available using the rjags method; to use these functions with other methods install (and specify using the module argument) the 'paretoprior' standalone module")	
			} 
		}
	}
	
	##  Get the data
	data <- checkdataformat(arg=data, block=character(0), auto=character(0), n.chains=NA, data.type=TRUE, evalscope=NULL)$combined
  
	##  Get the inits - evalscope is passed as a list with the character, so that it only a small reference is passed - it is only list.formatted if needed:
	inits <- checkdataformat(arg=inits, block=character(0), auto=character(0), n.chains=n.chains, data.type=FALSE, evalscope=list(data))
	n.chains <- inits$n.chains
	inits <- inits$combined
		
	data <- paste(data, "\n", sep="")
	if(gsub("[[:space:]]", "", data)=="")
		data <- ""
	if(identical(data, '') && runjags.getOption('nodata.warning'))
		warning('No data was specified or found in the model file so the simulation was run withut data', call.=FALSE)
	
	valid <- checkvalidforjags(data)	
	if(!valid$valid) stop(paste("The following problem was identified in the data provided:  ", valid$probstring, sep=""), call.=FALSE)		
	
	for(i in 1:length(inits)){
		valid <- checkvalidforjags(inits[i])	
		if(!valid$valid) stop(paste("The following problem was identified in the initial values provided for chain ", i, ":  ", valid$probstring, sep=""), call.=FALSE)	
	}
	
	
	
	if(length(grep('base::Mersenne-Twister', inits)>0) & as.numeric(jags.status$JAGS.major) < 2) warning('Using the RNG "base::Mersenne-Twister" (used by default for chain 4) may cause problems with restarting subsequent simulations using the end state of previous simulations due to a bug in JAGS version 1.x.  If you encounter the error "Invalid .RNG.state", please update JAGS to version 2.x and try again.  Or, you can change the random number generator by changing the .RNG.name to (for example) "base::Super-Duper" and remove the .RNG.state element of the list.', call.=FALSE)
	
	if(any(c("pd","full.pd","popt") %in% monitor) & (n.chains < 2 )) stop("The DIC, pD, full.pD and pOpt cannot be assessed with only 1 chain")
	if(any(c("pd","full.pd","popt","deviance") %in% monitor) & jags.status$JAGS.major < 2) stop('Support for the deviance, pD and popt monitors is no longer available for JAGS version 1.x.  Please update to JAGS version 3.x')
	
	
	# Combine model blocks, change \r to \n and get rid of double spacing:
	model <- paste(model, "\n", sep="", collapse="\n")
	model <- gsub("\r","\n",model)
	model <- gsub("\n\n","\n",model)
	model <- gsub("\n\n","\n",model)
	
	monitorcollapse <- ">\nmonitor set <"
	monitors <- paste("monitor set <", paste(monitor, collapse=monitorcollapse), ">\n", sep="")
	n.params <- length(monitor)
	params.names <- monitor
	
	
	if(!rjagsmethod){
		# Modules/factories and model setup/checking etc can be done by as.jags.runjags later
		temp.directory <- tempfile('runjagsdir')
		dir.create(temp.directory)
		cwd <- getwd()
		on.exit({
			setwd(cwd)
			unlink(temp.directory, recursive=TRUE)
			})
		setwd(temp.directory)

		cat(model, file="model.txt",sep="")  
		cat(data, file="data.txt",sep="")  
		for(i in 1:n.chains){
			cat(inits[i], file=paste("inits",i,".txt",sep=""),sep="")
		}	
		
		resetsyspath=resetbinpath <- FALSE
		if(.Platform$OS.type == "windows"){		
			currentsyspath <- Sys.getenv('PATH')
			if(!grepl(jags.status$libpaths$PATH,currentsyspath,fixed=TRUE)){
				Sys.setenv(PATH=paste(currentsyspath, ';', jags.status$libpaths$PATH, sep=''))
				resetsyspath <- TRUE
			}

			currentsysbinpath <- Sys.getenv('LTDL_LIBRARY_PATH')
			if(!grepl(jags.status$libpaths$LTDL_LIBRARY_PATH,currentsysbinpath,fixed=TRUE)){
				Sys.setenv(LTDL_LIBRARY_PATH=paste(currentsysbinpath, if(currentsysbinpath!='') ';', jags.status$libpaths$LTDL_LIBRARY_PATH, sep=''))
				resetbinpath <- TRUE
			}		
		}	
		
		scriptstring <- ""
		if(!identical(modules,"")) for(i in 1:length(modules)){
			if(modules[[i]][1]=="runjags") stop("The runjags module is only available using the rjags method; to use the functions provided with other methods install (and specify using the module argument) the 'paretoprior' standalone module")
			scriptstring <- paste(scriptstring, if(modules[[i]][2]=='FALSE') "un", "load ", modules[[i]][1], "\n", sep="")
		}
			
		if(!identical(factories,"")) for(i in 1:length(factories)){
			scriptstring <- paste(scriptstring, "set factory \"", factories[[i]][1], "\" ", if(factories[[i]][3]=='TRUE') "on" else "off", ", type(", factories[[i]][2], ")\n", sep="")
		}
	
		cat(scriptstring, "\nexit\n", file="script.cmd", sep="")  
		output <- system2(jags, stdout=TRUE, stderr=TRUE, stdin="script.cmd", wait=TRUE)
		if(grepl("file not found", paste(output,collapse="\n")) || grepl("error", tolower(deparse(paste(output,collapse="\n"))))){
			cat(output, sep="\n")
			stop("Error reading modules or factories (see output above for more details)")
		}
	
		scriptstring <- paste(scriptstring, "model in <\"model.txt\">\n", sep="")

		cat(scriptstring, "\nexit\n", file="script.cmd", sep="", append=FALSE)  
		output <- system2(jags, stdout=TRUE, stderr=TRUE, stdin="script.cmd", wait=TRUE)
		if(grepl("error", tolower(deparse(paste(output,collapse="\n"))))){
			cat(output, sep="\n")
			
			failedjagsmodel <- model
			class(failedjagsmodel) <- "runjagsmodel"
			assign("model", failedjagsmodel, envir=failedjags)
			stop("Error reading model (see output above for more details, and use failed.jags('model') to see model syntax with line numbers)")
		}
	
		if(data!=""){
			scriptstring <- paste(scriptstring, "data in <\"data.txt\">\n", sep="")
	
			cat(scriptstring, "\nexit\n", file="script.cmd", sep="", append=FALSE)  
			output <- system2(jags, stdout=TRUE, stderr=TRUE, stdin="script.cmd", wait=TRUE)
			if(grepl("error", tolower(deparse(paste(output,collapse="\n"))))){
				cat(output, sep="\n")
				class(data) <- 'runjagsdata'
				assign("data", data, envir=failedjags)
				
				stop("Error reading data (see output above for more details, and examine failed.jags('data') to see data file syntax with line numbers)")
			}
		}
	
		for(i in 1:n.chains){
			# Would have to compile the model to check inits, but we can check them as data:
			if(inits[i]!=""){
				cat(paste(scriptstring, "data in <\"inits",i,".txt\">\n", sep=""), "\nexit\n", file="script.cmd", sep="", append=FALSE)  
				output <- system2(jags, stdout=TRUE, stderr=TRUE, stdin="script.cmd", wait=TRUE)
				if(grepl("error", tolower(deparse(paste(output,collapse="\n"))))){
					cat(output, sep="\n")
					class(inits) <- 'runjagsinits'
					assign("inits", inits, envir=failedjags)
					
					stop(paste("Error reading initial values for chain ", i, " (see output above for more details, and use failed.jags('inits') to see init file syntax with line numbers)",sep=""), call.=FALSE)
				}
			}
		}
		
		if(resetsyspath) Sys.setenv(PATH=currentsyspath)
		if(resetbinpath) Sys.setenv(LTDL_LIBRARY_PATH=currentsysbinpath)
	}

	class(model) <- "runjagsmodel"
	class(inits) <- "runjagsinits"
	class(data) <- "runjagsdata"
	
	if(runjags.getOption('debug')>=10)
		cat('Monitors:  ', paste(monitor, collapse=', '), '\n', sep='')
	
	summary.pars <- getdefaultsummarypars()
	if(!identical(mutate, NA))
		summary.pars$mutate <- mutate
	# Otherwise the default value is NULL
	
	output <- list(mcmc=as.mcmc.list(lapply(1:n.chains, function(x) return(as.mcmc(NA)))), pd=NA, deviance.table=NA, deviance.sum=NA, end.state=inits, burnin=0, sample=0, thin=1, summary="", HPD="", hpd="", mcse="", psrf="", autocorr="", crosscorr="", dic="", trace=NA, density=NA, model=model, data=data, monitor=monitor, noread.monitor=character(0), modules=modules, factories=factories, response=response, fitted =fitted, residual=residual, method=method, method.options=getdefaultmethodoptions(), timetaken=0, summary.pars=summary.pars)
	class(output) <- "runjags"
	output <- checkvalidrunjagsobject(output)
	
	return(output)
}


setup.jagsfile <- function(model, n.chains=NA, data=NA, inits=NA, monitor=NA, modules=c(""), factories=c(""), mutate=NA, jags=runjags.getOption('jagspath'), method="simple", call.setup=TRUE, failincomplete=TRUE){
	
	# We may be passed some unevaluated function arguments so evaluate everything here:
	argnames <- names(formals(setup.jagsfile))
	for(i in 1:length(argnames)){
		success <- try(assign(argnames[i], eval(get(argnames[i]))), silent=TRUE)		
		if(inherits(success, 'try-error')){
			stop(paste("object '", strsplit(as.character(success),split="'",fixed=TRUE)[[1]][2], "' not found", sep=""), call.=FALSE)
		}
	}
	
	if(class(model)=="runjagsmodel")
		class(model) <- "character"
	if(identical(inits, list()))
		inits <- NA
	
	path <- model
	params <- read.jagsfile(path)

	autodata <- params$autodata
	autoinits <- params$autoinits
	maindata <- params$data
	maininits <- params$inits
	response <- params$response
	fitted <- params$fitted
	residual <- params$residual
	mainmutate <- params$mutate
	automutate <- params$automutate
	
	# Sort out mutate list, function or char:
	if(!identical(mutate,NA)){
		if((!is.na(mainmutate) || !is.na(automutate)) && runjags.getOption('blockignore.warning'))
			warning('Mutate functions specified in the model file or using #mutate# are ignored when a mutate argument is given', call.=FALSE)
	}else{
		if(!is.na(automutate) && !is.na(mainmutate)){
			warning('A mutate block was found in the model file, so #mutate# will be ignored', call.=FALSE)
			automutate <- NA
		}
		if(identical(mainmutate, NA)){
			mutate <- automutate
		}else{
			mutate <- mainmutate

			if(!identical(mutate, NA)){
				if(length(mutate)>1){
					warning("More than 1 mutate block was found in the file - the first was used and other(s) ignored", call.=FALSE)
					mutate <- mutate[1]			
				}
				s <- try(mutate <- eval(parse(text=mutate)))
				if(class(s)=='try-error'){
					warning('The mutate block provided could not be parsed and will be ignored', call.=FALSE)
					mutate <- NA
				}else{
					if(!class(mutate)%in%c('character','list','function')){					
						warning('The mutate block provided was not valid and will be ignored - it should either be a function, single character string or a list', call.=FALSE)
						mutate <- NA
					}
				}
			}
		}
	}	
	if(class(mutate)=='list'){
		if(class(mutate[[1]])=='character'){
			funtoget <- mutate[[1]]
		}else{			
			if(!is.function(mutate[[1]]))
				stop('The fist element of the list given as the mutate argument must either be a function or character matching a function', call.=FALSE)
			funtoget <- mutate[[1]]
		}
	}else{
		funtoget <- mutate
	}
	# Try finding the function now rather than later (but don't save it):
	if(!identical(funtoget, NA) && !is.function(funtoget)){
		if(!is.character(funtoget) || length(funtoget)!=1)
			stop('The mutate argument must either be a function, character string matching a function, or a list with the first element a function or character string', call.=FALSE)
		s <- try(funfound <- get(funtoget))
		if(class(s)=='try-error' || !is.function(funfound))
			stop('The character string argument supplied to mutate did not match a valid function', call.=FALSE)
	}
	
	
	oldmodules <- as.list(checkmodfact(modules, 'module'))
  	newmodules <- as.list(checkmodfact(params$modules, 'module'))
	modules <- checkmodfact(c(oldmodules,newmodules),'module')

	oldfactories <- as.list(checkmodfact(factories, 'factory'))
  	newfactories <- as.list(checkmodfact(params$factories, 'factory'))
	factories <- checkmodfact(c(oldfactories,newfactories),'factory')
	
	if(params$model=="model{\n\n}\n") stop("No valid model was specified or found in the model block")
	outmodel <- params$model
	
	##  Get the data
	if(inherits(data, c('runjagsdata','character'))){
		if((!identical(maindata, NA) || !identical(autodata, NA)) && runjags.getOption('blockignore.warning'))
			warning('Data specified in the model file or using #data# are ignored when a character string is given as the argument to data', call.=FALSE)
		maindata <- NA
		autodata <- NA
		outdata <- data
	}else{
		outdata <- checkdataformat(arg=data, block=maindata, auto=autodata, n.chains=NA, data.type=TRUE, evalscope=NULL)$combined
	}
	
	if(all(sapply(inits,class)%in%c('runjagsinits','character'))){
		if((!identical(maininits, NA) || !identical(autoinits, NA)) && runjags.getOption('blockignore.warning'))
			warning('Inits specified in the model file or using #inits# are ignored when a character string is given as the argument to inits', call.=FALSE)
		maininits <- NA
		autoinits <- NA
		outinits <- unlist(inits)
		n.chains <- length(inits)
	}else{
		##  Get the inits - evalscope is passed as a list with the character, so that it only a small reference is passed - it is only list.formatted if needed:
		outinits <- checkdataformat(arg=inits, block=maininits, auto=autoinits, n.chains=n.chains, data.type=FALSE, evalscope=list(outdata))
		n.chains <- outinits$n.chains
		outinits <- outinits$combined
	}
	
	if(outdata==''){
		outdata <- NA
#		if(failincomplete && runjags.getOption('blockignore.warning'))
#			warning("The model was run without data since no data was provided or found in the model block", call.=FALSE)
# I can deal with a lack of data better by looking at the samplers
	}
	
	if(!identical(monitor, NA) && !identical(params$monitor, NA)){
		if(runjags.getOption('blockignore.warning'))
			warning('Monitors specified in the model file are ignored when a monitors argument is given', call.=FALSE)
		params$monitor <- NA
	}
	monitor <- (unique(na.omit(c(monitor, params$monitor))))
	outmonitor <- monitor[monitor!='']
	if(length(outmonitor)==0 && failincomplete)
		stop("No monitors were specified or found in the model block", call.=FALSE)
	
	# Ignored pretty much everywhere - the MCMC chains are returned alphabetically ALWAYS - NOT TRUE!!!!!!
#	outmonitor <- sort(outmonitor)
#	if(any(outmonitor=='resid.sum.sq'))
#		outmonitor <- c(outmonitor[outmonitor!='resid.sum.sq'], 'resid.sum.sq')
#	if(any(outmonitor=='deviance'))
#		outmonitor <- c(outmonitor[outmonitor!='deviance'], 'deviance')
	
	lengths <- lapply(params, length)
	if(any(lengths==0) && failincomplete)
		stop(paste("No ", paste(names(lengths[lengths==0]), collapse=" or "), " blocks or tags were found", sep=""))
	
	if(runjags.getOption('debug')>=10){
		if(identical(outdata, NA))
			cat('No data\n')
		else
			cat(outdata)
		cat(outinits)
	}
  
	if(call.setup){
		return(setup.jags(model=outmodel, monitor = outmonitor, data=outdata, n.chains=n.chains, inits = outinits, modules=modules, factories=factories, response=response, fitted =fitted, residual=residual, jags = jags, method=method, mutate=mutate))
	}else{
		return(list(model=outmodel, monitor = outmonitor, data=outdata, n.chains=n.chains, inits = outinits, modules=modules, factories=factories, response=response, fitted =fitted, residual=residual, jags = jags, mutate=mutate))
	}
	
}


checkdataformat <- function(arg, block, auto, n.chains=NA, data.type=TRUE, evalscope=NULL){
	
	# This function will either return a character (if given a single character), or a list (if other combinations)
	# This all works with inits too - need to specify evalscope (data) and n.chains though	
	
	# Used a couple of times - character has to be included for the benefit of .RNG.name:
	permitted.formats <- c('factor','numeric','integer','double','logical','character','array','matrix')
	
	
	# Otherwise unevalauted functions are a problem:
	arg <- eval(arg)
	block <- eval(block)
	auto <- eval(auto)
	
	if(any(is.na(auto)))
		auto <- na.omit(auto)
	
	if(data.type) string <- 'data' else string <- 'inits'
	
	if(identical(arg, NA) || identical(arg, '') || is.null(arg)){
		arg <- character(0)
	}
	if(identical(block, NA) || identical(block, '') || is.null(block)){
		block <- character(0)
	}
	if(identical(auto, NA) || identical(auto, '') || is.null(auto)){
		auto <- character(0)
	}
	
	stopifnot(length(n.chains)==1)
	
	# Set n.chains:
	if(data.type){	  
    	n.chains <- 1
	}else{
		# If working with inits, try to guess n.chain using:
		# 1) The length of lists of lists (or length of chars or list of chars)
		# 2) The number of blocks
		# If they conflict, take the max and give a warning

		if(is.na(n.chains)){
			
      	  n.chain.block=n.chain.list <- NA
			if(!length(block)==0){
				n.chain.block <- max(length(block),1)
			}					
			if(inherits(arg, c('character', 'runjagsinits', 'list')) && length(arg)>0 && all(sapply(arg, inherits, what=c('character','list','data.frame','runjagsinits','environment')))){
				n.chain.list <- max(length(arg),1)
			}
			# If they are both specified and don't match lengths:
			if(!is.na(n.chain.block) && !is.na(n.chain.list) && n.chain.block!=n.chain.list && runjags.getOption('inits.warning')){
				warning("The number of initial value blocks specified in the model file did not match the length of the inits argument - using the larger number of chains", call.=FALSE)
			}
      # Horrible hack to make n.chains -10 if not already set:
			n.chains <- max(n.chain.block, n.chain.list, -10, na.rm=TRUE)
			stopifnot(!is.na(n.chains))
      stopifnot(n.chains==-10 || n.chains >=0)			
			if(n.chains==-10){
				n.chains <- 2
				if(runjags.getOption('inits.warning'))
					warning("No initial value blocks found and n.chains not specified: 2 chains were used", call.=FALSE)
			}
		}
	}
	
  	stopifnot(length(n.chains)==1 && n.chains>0)
  	
	# First deal with a funcion - may need to assign the evalscope if it is inits:
	if(inherits(arg, 'function')){
		if(data.type){
			if(!is.null(formals(arg)))
				stop(paste("The function provided for ", string, " must take exactly zero arguments", sep=""), call.=FALSE)
			s <- try(newarg <- arg())
			if(class(s)=="try-error")
				stop(paste("The following error was encountered using the function supplied for ", string, ":  ", s[1], sep=''), call.=FALSE)
			arg <- list(newarg)
		}else{
			# Need to have been passed a list with a single reference to the character (too expensive to pass a massive char)
			stopifnot(class(evalscope)=='list' && length(evalscope)==1)
			evalscope <- list.format(evalscope[[1]])
			
			### In case the same function environment is being used by many simulations (this is definitely necessary):
			# Copy the function:
			localinits <- arg
			# Create a new environment with the parent as the environment of the original function:
			newenv <- new.env()
			parent.env(newenv) <- environment(arg)
			# This ensures the following data is first on the search path but anything else around from this environment is also available:
			assign('data', evalscope, newenv)
			# Set the environment to the new local environment:
			environment(localinits) <- newenv
		  
			if(!length(formals(localinits))%in%c(0,1))
			stop("The function provided for 'inits' must take either zero arguments or one argument representing the chain number",call.=FALSE)
      
      	  	s <- try({
			if(is.null(formals(localinits))){
				arg <- lapply(1:n.chains, function(x) return(localinits()))
			}else{			
				arg <- lapply(1:n.chains, function(x) return(localinits(x)))			
			}
			})
			if(class(s)=="try-error") stop(paste("The following error was encountered using the function supplied for 'inits':  ", s[1], sep=''), call.=FALSE)
		}
	}
	
	# If working with data need to make this into a length 1 list:
	if(data.type){
		arg <- list(arg)
	}
  
	# Now convert single specified arguments into a list:
	if(inherits(arg, c('runjagsdata', 'runjagsinits', 'character')))
		arg <- as.list(arg)
	if(inherits(arg, c('data.frame', 'environment')))
		arg <- list(arg)
  
  	# Now catch a single named list rather than a list of lists:
	if(inherits(arg,'list') && !is.null(names(arg)) && all(sapply(arg, class)%in%permitted.formats))
		arg <- list(arg)
	
	# If arg doesn't match n.chains recycle it:
	if(length(arg)>0 && length(arg)!=n.chains){
		stopifnot(!data.type)
		if(runjags.getOption('inits.warning') && !all(sapply(arg,class)=='environment'))
			warning("The length of the initial values argument supplied found does not correspond to the number of chains specified.  Some initial values were recycled or ignored.", call.=FALSE)
	
		usechain <- rep(1:length(arg), length.out=n.chains)
		arg <- lapply(usechain, function(x){
      return(arg[[x]])
      })
	}
	# Otherwise make it an emptylist if it is length 0:
	if(length(arg)==0)
		arg <- replicate(n.chains,NULL)
	
	stopifnot(class(block)=='character')
	block <- as.list(block)
	# If block doesn't match n.chains recycle it:
	if(length(block)>0 && length(block)!=n.chains){
		stopifnot(!data.type)
		if(runjags.getOption('inits.warning'))
			warning("The number of initial value blocks found does not correspond to the number of chains specified.  Some initial values were recycled or ignored.", call.=FALSE)
	
		usechain <- rep(1:length(block), length.out=n.chains)
		block <- lapply(usechain, function(x) return(block[[x]]))
	}
	
	# Now we are guaranteed to have a list of things to work with for arg (list may be of length 0)
	# But we need to make sure each component of the list is viable (length arg should be length chains):
	for(i in 1:length(arg)){
		if(inherits(arg[[i]], 'data.frame'))
			arg[[i]] <- as.list(arg[[i]])
		if(inherits(arg[[i]], c("runjagsdata","runjagsinits")))
			class(arg[[i]]) <- "character"
		if(!is.null(arg[[i]]) && !(inherits(arg[[i]], c('character','list','environment'))))
			stop(paste('Unrecognised format for the argument specified as ', string, if(!data.type) paste(' (chain ', i, ')', sep=''), ' - it must either be in R dump format (see dump.format(), or a named list, data frame or environment, or a function returning one of these types', sep=''), call.=FALSE)

		if(inherits(arg[[i]], 'list') && (is.null(names(arg[[i]])) || any(names(arg[[i]])=='') || length(unique(names(arg[[i]])))!=length(arg[[i]])))
			stop(paste('The list specified as the argument for ', string, if(!data.type) paste(' (chain ', i, ')', sep=''), ' was not a fully named list', sep=''), call.=FALSE)		
	}
	stopifnot(length(arg)==n.chains)
	
	# Loop through chains to find auto, and any naming conflicts:
	globnames <- ls(envir=.GlobalEnv, all.names=TRUE)
	toreturn <- character(n.chains)
	for(c in 1:n.chains){
		
		# If the function has evalauted to a character, remove the auto founds:
		if(is.character(arg[[c]]))
			auto <- character(0)

		# If we have any auto find the auto grabbed variables:
		autofound <- auto
		autonames <- character(0)
		if(length(auto)>0){
			env <- arg[[c]]
			if(inherits(env, c('list', 'data.frame')))
				env <- as.environment(as.list(env))
			if(class(env)!='environment')
				env <- .GlobalEnv
			
			envnames <- ls(envir=env, all.names=TRUE)
			geninits <- vector('list', length=length(auto))
			
			for(i in 1:length(auto)){	
				feedtxt <- paste(string, " specified", sep="")
				chaintxt <- if(!data.type) paste(' (for chain ', i, ')', sep='') else ''
				if(any(envnames==auto[i])){
					temp <- get(auto[i], envir=env, inherits=TRUE)
				}else{
					if(any(globnames==auto[i])){
						temp <- get(auto[i], envir=.GlobalEnv, inherits=TRUE)
						feedtxt <- 'in the global environment'
						chaintxt <- ''
					}else{
						stop(paste(auto[i], " was not found in the ", feedtxt, " or the global environment", sep=""), call.=FALSE) 
					}
				}
				
				if(class(temp)=="function"){
					success <- suppressWarnings(try(temp <- temp(c), silent=TRUE))
					if(inherits(success, 'try-error') && !data.type){
					  success <- suppressWarnings(try(temp <- temp(), silent=TRUE))
					}
					if(inherits(success, 'try-error'))					
            			stop(paste('A function matching the variable name "', auto[i], '" was found in the ', feedtxt, ' but cannot be executed with 0', if(!data.type) ' or 1', ' arguments', sep=''), call.=FALSE)
					
					if(inherits(temp, c('data.frame','environment')))
						temp <- as.list(temp)
					
					if(data.type && class(temp)=='list')
						stop(paste('A function matching the variable name "', auto[i], '" was found in the ', feedtxt, ' but returns a list which is not allowed for data - it should return either a factor, numeric, integer or logical', sep=''), call.=FALSE)											
					
					if(is.null(temp))
						stop(paste('A function matching the variable name "', auto[i], '" was found in the ', feedtxt, ' but returns NULL', sep=''), call.=FALSE)
					
					if(!class(temp) %in% c(permitted.formats, 'list') || (class(temp)=='list' && !all(sapply(temp,class)%in%permitted.formats)))
            			stop(paste('A function matching the variable name "', auto[i], '" was found in the ', feedtxt, ' but returns an object of class ', class(temp), ' when executed (this should be a factor, numeric, integer, logical, or a list of one of these types)', sep=''), call.=FALSE)
				}
				
				if(class(temp)=="list"){
					if(data.type)
						stop(paste('The variable name "', auto[i], '" was found in the ', feedtxt, ' but it is a list which is not allowed for data - it should be a factor, numeric, integer or logical', sep=''), call.=FALSE)					
					
					# Allow recycling values if the list is of length 1 (this is equivalent to a non-list anyway):
					if(length(temp) < c){
						if(length(temp)==1)
							temp <- temp[[1]]
						else
							stop(paste('A list matching the variable name "', auto[i], '" was found in the ', feedtxt, ' but has a different length to the number of chains specified', sep=''), call.=FALSE)
					}else{
						temp <- temp[[c]]
					}
					if(!class(temp) %in% permitted.formats)
	        			stop(paste('A list matching the variable name "', auto[i], '" was found in the ', feedtxt, ' but the element corresponding to chain ', c, ' is of class ', class(temp), ' - it should be a factor, numeric, integer or logical', sep=''), call.=FALSE)
					
					if(is.null(temp))
						stop(paste('A list matching the variable name "', auto[i], '" was found in the ', feedtxt, ' but the element corresponding to chain ', c, ' is NULL', sep=''), call.=FALSE)
				}
				
				# necessary to remove compound listing somehow introduced by initlist function or something:
				while(class(temp)=="list")
					temp <- temp[[1]]
				
				if(!class(temp) %in% permitted.formats)
        			stop(paste('An object matching the variable name "', auto[i], '" was found in the ', feedtxt, chaintxt, ' but is of class ', class(temp), ' - it should be a factor, numeric, integer or logical', sep=''), call.=FALSE)
				
				geninits[[i]] <- temp

			}
			names(geninits) <- auto			
			autofound <- dump.format(geninits, checkvalid=TRUE)
			
		}else{
			
      	  	# Otherwise if there aren't any auto then just take the whole list or data frame, but ignore an env or character:
			if(inherits(arg[[c]], c('list','data.frame'))){
				autofound <- dump.format(arg[[c]], checkvalid=TRUE)
				autonames <- names(arg[[c]])
			}else{
				# This will happen for an env or character function return:
				autofound <- ''
			}
			# Note that character specified arg goes to block UNLESS it is the result of a function!!!
		}
		
		blockfound <- if(length(block)==0) '' else block[[c]]
		if(length(blockfound)==0)
			blockfound <- ''
		
		argfound <- ''
		if(is.character(arg[[c]])){
			argfound <- arg[[c]]
		}
		
		# If we are using 2 or 3 different ones:
		if(sum((argfound!='') + (blockfound!='') + (autofound!=''))>1){
			blocknames <- names(list.format(blockfound, checkvalid=TRUE))
			argnames <- names(list.format(argfound, checkvalid=TRUE))
			
			allnames <- c(argnames, blocknames, autonames)
			fnames <- factor(allnames)
			occurances <- tabulate(fnames)
			if(any(occurances)>1){
				problems <- levels(fnames)[which(occurances>1)]
				stop(paste('The following variables have a naming conflict for different sources of ', string, if(!data.type) paste(' (chain ', c, ')', sep=''), ': ', paste(problems, collapse=', '), sep=''), call.=FALSE)
			}
		}
		
		toreturn[c] <- paste(autofound, blockfound, argfound, '\n', sep='\n', collapse='\n')
		
		if(identical(toreturn[c], NA))
			toreturn[c] <- ""
	
		# Change \r to \n and get rid of double spacing and leading \n:
		toreturn[c] <- gsub("\r","\n",toreturn[c])
		while(grepl("\n\n", toreturn[c]))
			toreturn[c] <- gsub("\n\n","\n",toreturn[c])
		toreturn[c] <- gsub("^[[:space:]]*\n", "", toreturn[c])
		
		if(gsub("[[:space:]]", "", toreturn[c])=="")
			toreturn[c] <- ""
		
	}
		
	return(list(combined=toreturn, n.chains=n.chains))	
}


setup.JAGS <- setup.jags