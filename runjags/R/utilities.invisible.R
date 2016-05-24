getRvers <- function(){
	vers <- R.version$version.string
	version.string <- gsub("\\(.*?\\)","",vers,perl=FALSE)
	return(getvers(version.string))
}

getvers <- function(version.string){
	vers <- gsub("[[:space:]]","",gsub("[[:alpha:]]","",version.string))
	vers <- gsub("[[:punct:]]", "", gsub("."," ",vers, fixed=TRUE))
	vers <- strsplit(vers," ",fixed=TRUE)[[1]]
	version <- (10^((length(vers)-1):0))^3 * as.numeric(vers)
	version <- sum(version)
	return(version)
}

getjagsnames <- function(targets){
	
	retval <- unlist(lapply(1:length(targets), function(i){
	
		x <- targets[[i]]
		n <- names(targets)[i]
		new <- as.list(x)

		if(length(x)==1){
			names(new) <- n
			return(new)
		}
		if(!is.array(x)){
			names(new) <- paste(n, "[",1:length(x),"]",sep="")
			return(new)
		}else{
			dims <- apply(expand.grid(lapply(dim(x), function(x) return(1:x))),1,paste,collapse=",")
			names(new) <- paste(n,"[",dims,"]",sep="")
			return(new)
		}
		stop(paste("Unsupported argument type:",class(x)))
	}))
	
	return(retval)
}

# Currently used only by glm.template - rjags and JAGS now respect initial ordering of variables:
alphabeticalvars <- function(x, always.last=c('resid.sum.sq','deviance')){
	
	stopifnot(class(x)=='character')
	stopifnot(class(always.last)=='character')
	stopifnot(all(!is.na(always.last)))

	splits <- strsplit(gsub(']','',x,fixed=TRUE),'[',fixed=TRUE)	
	mnames <- sapply(splits, function(x) return(x[[1]]))

	indices <- lapply(splits, function(x){
		if(length(x)==1)
			return(numeric(0))
		else
			return(as.numeric(strsplit(x[2],',')[[1]]))
		})

	lengths <- sapply(indices,length)
	maxind <- max(lengths)+1

	sortmat <- matrix(0, nrow=length(x), ncol=maxind)
	for(i in 1:length(x)){
		if(lengths[i]>0)
			sortmat[i,2:(lengths[i]+1)] <- indices[[i]]
	}
	sortmat[,1] <- rank(mnames, ties.method='average')
	
	# If any match the always lasts chars:
	if(length(always.last)>0){
		for(i in 1:length(always.last)){
			sortmat[which(x==always.last[i]),1] <- nrow(sortmat)+i
		}
	}
	
	arglist <- lapply(1:ncol(sortmat), function(x) return(as.numeric(sortmat[,x])))
	return(x[do.call('order', args=arglist)])
	
}

getarraynames <- function(targets){
	
	vars <- gsub('\\[[[:print:]]*\\]','',names(targets))
	varnames <- unique(vars)
	indexes <- gsub("[[:alpha:][:punct:]]", "", gsub(","," ",names(targets)))
	indexes[indexes==""] <- 1
	dimensions <- lapply(indexes, function(y){		
		if(!grepl(' ',y)) return(as.numeric(y)) else return(as.numeric(strsplit(y, split=" ",fixed=TRUE)[[1]]))		})
	
	if(any(sapply(varnames, function(x) return(grepl(',',x)))))
		stop("Invalid symbol in variable name")
	
	retlist <- lapply(varnames, function(x){
		ds <- sapply(dimensions[vars==x],length)
		stopifnot(all(ds==ds[1]))
		thisdim <- sapply(dimensions[vars==x], function(y) return(y))
		if(is.null(dim(thisdim))) dim(thisdim) <- c(1,length(thisdim))
    dims <- apply(thisdim,1,max)
		newarr <- targets[vars==x]
		dim(newarr) <- dims
		return(newarr)
	})
	names(retlist) <- varnames
	
	return(retlist)
	
}

makerunjagsobject <- function(combinedoutput, summarise, summaryargs, burnin, sample, thin, model, data, monitor, noread.monitor, modules, factories, response, residual, fitted, method, method.options, timetaken, silent=FALSE){
	
	starttime <- Sys.time()
	
	summaryargs <- getargs('add.summary', summaryargs)
	evalsumargs <- lapply(summaryargs, eval, envir = parent.frame(), enclos=baseenv())   # some get evaluated to NULL so must set up a new vector to avoid indexing problems - this scoping means check the parent frame first, but find runjags.getOption in runjags environment
	
	combinedoutput <- combinedoutput[names(combinedoutput)%in%c('mcmc','deviance.table','deviance.sum','pd','end.state','samplers')]
	stopifnot(identical(names(combinedoutput), c('mcmc','deviance.table','deviance.sum','pd','end.state','samplers')))
	
	combinedoutput <- c(combinedoutput, list(burnin=burnin, sample=niter(combinedoutput$mcmc), thin=thin, model=model, data=data, monitor=monitor, noread.monitor=noread.monitor, modules=modules, factories=factories, response=response, residual=residual, fitted=fitted, method=method, method.options=method.options, timetaken=timetaken, runjags.version=c(runjagsprivate$runjagsversion, R.Version()$version.string, .Platform$OS.type, .Platform$GUI, .Platform$pkgType, format(Sys.time()))))
	class(combinedoutput) <- 'runjags'
	combinedoutput$summary.pars <- getdefaultsummarypars()
	
	if(is.na(summarise)) summarise <- TRUE
	# If too many vars (>50?), don't summarise or plot unless forced to:
	if(summarise && !runjags.getOption('force.summary') && nvar(combinedoutput$mcmc)>50){
		# Don't produce the warning - it breaks one of the tests (fb('Check that summary statistics arent calculated for >50 variables (and it can be overridden)') for a start
#		if(!silent)
#			warning("Summary statistics were not produced as there are >50 monitored variables - to override this behaviour see ?add.summary and ?runjags.options", call.=FALSE)
		summarise <- FALSE
	}
		
	# Call function to calculate summary statistics and plots etc:
	if(summarise){
		success <- try({
		summaries <- do.call('add.summary', c(list(runjags.object=combinedoutput), summaryargs))
		}, silent=TRUE)
		if(!silent && inherits(success, 'try-error')){
			warning(paste('The following error occured while calculating summary statistics:\n  ', gsub('\n$','',as.character(success)), sep=''),call.=FALSE)
			message <- "An error occured while calculating summary statistics - see ?add.summary"
			summaries <- list(summary=message, HPD=message, hpd=message, mcse=message, psrf=message, autocorr=message, crosscorr=message, stochastic=message, trace=message, density=message, hist=message, ecdfplot=message, key=message, acplot=message, ccplot=message, summaries=NULL, summary.available=FALSE, summary.pars=evalsumargs)
		}
	}else{
		message <- "Summary statistics are not stored internally when summarise=FALSE - see ?add.summary"
		summaries <- list(summary=message, HPD=message, hpd=message, mcse=message, psrf=message, autocorr=message, crosscorr=message, stochastic=message, trace=message, density=message, hist=message, ecdfplot=message, key=message, acplot=message, ccplot=message, summaries=NULL, summary.available=FALSE, summary.pars=evalsumargs)
	}
	if(any(c('dic','ped') %in% combinedoutput$monitor)){
		summaries <- c(summaries, list(dic=runjags.dic(deviance.table=combinedoutput$deviance.table, deviance.sum=combinedoutput$deviance.sum, mcmclist=combinedoutput$mcmc)))
	}else{
		summaries <- c(summaries, list(dic="DIC statistics have not been stored internally - see ?extract.runjags"))
	}

	combinedoutput <- combinedoutput[! names(combinedoutput)%in%names(summaries)]
	combinedoutput <- c(combinedoutput, summaries)
	class(combinedoutput) <- 'runjags'
	
	class(combinedoutput$model) <- 'runjagsmodel'
	class(combinedoutput$data) <- 'runjagsdata'
	class(combinedoutput$end.state) <- 'runjagsinits'

	class(combinedoutput) <- 'runjags'
	
	# Check to see if we have identical RNG states after the run:
	n.chains <- nchain(combinedoutput$mcmc)
	if(!silent && n.chains > 1 && any(combinedoutput$end.state[2:n.chains]==combinedoutput$end.state[1])){
		warning("Identical RNG states have been produced for multiple chains - this probably means the chains are identical.  Try again using different parameter values and/or RNG samplers as inits.", call.=FALSE)
	}
	
	# If the model was compiled previously with lecuyer loaded but not listed as a module the lecuyer module will still not be listed, so add it:
	if(any(grepl('lecuyer::RngStream', combinedoutput$end.state, fixed=TRUE)))
	  combinedoutput$modules <- c(combinedoutput$modules, list(c("lecuyer","TRUE")))
	
	combinedoutput$modules <- checkmodfact(combinedoutput$modules,'module')
	combinedoutput$factories <- checkmodfact(combinedoutput$factories,'factory')
	
	# Add the summarise time to the time taken:
	combinedoutput$timetaken <- combinedoutput$timetaken + difftime(Sys.time(), starttime)
	
	return(combinedoutput)
}

swcat <- function(...){
	
	if(!runjags.getOption('silent.runjags')){
		pargs <- list(...)
		pasted <- do.call(paste, pargs)
		pasted <- gsub('\r', '\n', pasted)
	
		# White space is destroyed by strwrap so preserve \n by splitting on them (and append a ' ' [which is removed by strwrap anyway] to preserve any trailing \n)
		pasted <- unlist(strsplit(paste(pasted,' ',sep=''), '\n'))
		pasted <- strwrap(pasted)
		cat(paste(pasted, collapse='\n'))
	}
}

versionmatch <- function(required, actual){
	actual <- as.character(actual)
	
	matched <- FALSE
	for(r in as.character(required)){
		# Default:
		type <- "gteq"
		# If only an equals match precise version:
		if(grepl("=", r, fixed=TRUE)) type <- "eq"
		# Greater than takes precedence:
		if(grepl(">", r, fixed=TRUE)) type <- "gt"
		# Greater than or equal also possible:
		if(grepl(">=", r, fixed=TRUE)) type <- "gteq"	
		r <- gsub(">|=", "", r)
		if((compareVersion(actual, r)==0) & (type=="eq" | type=="gteq")){
			matched <- TRUE
		}
		if((compareVersion(actual, r)==1) & (type=="gt" | type=="gteq")){
			matched <- TRUE
		}
	}
	
	return(matched)	
}


winbugs.extract.big <- function(find, string, remove.list=TRUE){

	stopifnot(length(find)==1)
	stopifnot(length(string)==1)
	# There will not be any # except for the bugsdata jagsdata etc comments with the rest of the line cleared

	if(identical(string, '\n') || identical(string, '\r') || identical(string, '') || !grepl(find, string, ignore.case=TRUE))
	  return(list('', listfound=FALSE))
  
	# Make sure \r are now \n and \t are :
	string <- gsub('\r', '\n', string)
	string <- gsub('\t', ' ', string)

	# Remove white space at the beginning of lines but make sure a \n starts the string otherwise we are in trouble:
	string <- gsub('\n[[:space:]]*', '\n', paste('\n', string, sep=''))


	#########   Step 1:  separate blocks on the find string:

	# Now for the find to be valid it must be following a \n and before a {
	if(!grepl(paste('\n', find, '[[:space:]]*\\{',sep=''), string, ignore.case=TRUE))
	  return(list('', listfound=FALSE))

	# Now we have to turn these into lowercase because strsplit doesnt have an ignore.case argument:
	string <- gsub(paste('\n', find, '[[:space:]]*\\{',sep=''), paste('\n', find, ' \\{',sep=''), string, ignore.case=TRUE)

	# Now we can split on the find as we know there is at least 1, so 2 fragments after splitting:
	sepstrings <- strsplit(string, paste('\n', find, ' \\{',sep=''))[[1]]
	stopifnot(length(sepstrings)>=2)
	sepstrings <- sapply(sepstrings, function(x) return(paste(find, ' {\n', x, sep='')))[2:length(sepstrings)]
	names(sepstrings) <- NULL

	# Now the start of the strings is correct, but there may be crap after the closing } (but we might have a loop inside the string to keep):
	cutstrings <- lapply(sepstrings, function(s){
		# Use gregexpr to find the locations of all { and }:
		openbraces <- as.numeric(gregexpr('\\{', s)[[1]])
		closebraces <- as.numeric(gregexpr('\\}', s)[[1]])
		# Need to not count { } on commented lines, but also leave comment lines as they are - this finds where commented { } are:
		commented <- gregexpr("#(.*?)\\{", s, perl=TRUE)[[1]]
		openbraces <- openbraces[!openbraces %in% (as.numeric(commented)+(attr(commented, 'match.length')-1))]
		commented <- gregexpr("#(.*?)\\}", s, perl=TRUE)[[1]]
		closebraces <- closebraces[!closebraces %in% (as.numeric(commented)+(attr(commented, 'match.length')-1))]
		
		# Sanity check - an open brace should be right at the start:
		stopifnot(min(openbraces) < min(closebraces))
		# If we have different numbers of close brackets than open brackets, bail out:
		if(length(openbraces) > length(closebraces))
			stop(paste('There was an unmatched number of { and } found following a ', find, ' block - ensure that the text "', find, '{" is not used inside another block', sep=''), call.=FALSE)
		# This will tell us which closebrace corresponds to the last open brace:
		bracesdone <- which(sapply(1:length(closebraces), function(b) return(sum(openbraces < closebraces[b])==b)))[1]
		# And the matching string to return (from after the opening to before the closing bracket):
		return(substr(s, openbraces[1]+1, closebraces[bracesdone]-1))
	})
	cutstrings <- sapply(cutstrings, paste, '\n', sep='')
	# Remove leading whitespace again:
	cutstrings <- sub('^[[:space:]]*','',cutstrings)
	# And trailing whitespace again:
	cutstrings <- sub('[[:space:]]*$','',cutstrings)
	
	listsfound <- grepl('[[:space:]]*list[[:space:]]*\\(', cutstrings, ignore.case=TRUE)

	return(list(cutstrings, listfound=listsfound))

	#### We no longer separate/remove list here but the code was:
	
	########  Step 2:  separate and remove list( )

	# # If we don't have any or aren't removing them just bail:
	# listsfound <- grepl('[[:space:]]*list[[:space:]]*\\(', cutstrings, ignore.case=TRUE)
	# if(!any(listsfound) || !remove.list)
	#   return(list(cutstrings, listfound=listsfound))
	# 
	# # Otherwise we need to remove the list (
	# cutstrings <- gsub(paste('[[:space:]]*list[[:space:]]*\\(',sep=''), ' ', cutstrings, ignore.case=TRUE)
	# # And the last closebracket found (assume this is correct):
	# tozap <- cutstrings[listsfound]
	# tozap <- sapply(tozap, function(s){
	#   closebraces <- as.numeric(gregexpr('\\)', s)[[1]])
	#   toremove <- closebraces[length(closebraces)]
	#   substr(s, toremove, toremove) <- ' '
	#   return(s)
	# })
	# cutstrings[listsfound] <- tozap
	# 
	# # Remove leading whitespace again:
	# cutstrings <- sub('^[[:space:]]*','',cutstrings)
	# # And trailing whitespace again:
	# cutstrings <- sub('[[:space:]]*$','',cutstrings)
	# 
	# # And return cutstrings:
	# return(list(cutstrings, listfound=listsfound))

}


sortjagsvsbugs <- function(maindata, data.type=TRUE){
	
	# This is all worded for data, but it applies equally to inits as well
	feedtxt <- if(data.type) 'A data' else 'An inits'
	
	# Note that the list( ) will already have been removed - datalistfound tells us if it was there:	
	datalistfound <- maindata[[2]]
	maindata <- maindata[[1]]
	
	transwarngiven <- FALSE
	
	modelret <- ''
	for(i in 1:length(maindata)){
		
		if(maindata[i]=='')
			next
		
		# Manual override mechanism for data types
    	tofinds <- tolower(maindata[i])
		bugsdata <- grepl("#bugsdata#", tofinds)   	# =fixed
		jagsdata <- grepl("#rdata#", tofinds)  	# =fixed
		modeldata <- grepl("#modeldata#", tofinds)
		
		# Removes any  character between '#' and a new line
		maindata[i] <- gsub("#(.*?)\n", "\n", maindata[i], perl=TRUE)
		
		# Fixed means it is to be read as data (jags or bugs column/row major order), model means part of the model:
		# models don't contain .Dim, c(, ' or " (I think)
		likelyfixed <- any(c(grepl(".Dim", maindata[i]), datalistfound[i], grepl("c\\(", maindata[i]), grepl("'", maindata[i]), grepl('"', maindata[i])))
		# fixed data doesn't contain [ or for(
		likelymodel <- data.type && any(c(grepl("\\[", maindata[i]), grepl("for *\\(", maindata[i])))
		# Since I'm removing comments I think it's safe to assume anything with = is fixeddata?, anything with .Dim is fixeddata, anything with list( is fixeddata, anything with c( is fixeddata, anything with [ is modeldata, anything with for( is modeldata.  Check all these and if we get a disparity then a warning?
    # Manual overrides are done above		
		
		if(!data.type && modeldata){
			warning('Ignoring an invalid #modeldata# specification in the initial values', call.=FALSE)
			modeldata <- FALSE
		}
		
		# If we don't have an override:
		if(!any(bugsdata, jagsdata, modeldata)){
			if((likelyfixed && likelymodel) || (!likelyfixed && !likelymodel)){
				if(runjags.getOption('debug'))
					swcat(feedtxt, ' type was not recognised - assuming fixed type ', if(datalistfound[i]) '(BUGS format)' else '(R format)', '\n', sep='')

				#print("try parsing here?  Have an option to default data type if unkown?")
				#for(i in 1:length(maindata)){
				#	err <- try(eval(parse(text=maindata[i])))
				#}

				# For now assume that default is jags type data:
				modeldata <- FALSE
				bugsdata <- datalistfound[i]
				jagsdata <- !datalistfound[i]
			}else{
				modeldata <- likelymodel
				bugsdata <- !likelymodel && datalistfound[i]
				jagsdata <- !likelymodel && !datalistfound[i]
			}
		}
		
		if(modeldata){
			if(modelret!='')
				stop('Multiple model style data blocks were found in the model - this is not allowed', call.=FALSE)
			# Paste this data block onto the model and remove it from the data:
			modelret <- paste('data{\n', maindata[i], '\n}\n\n', sep='')
			maindata[i] <- ''
		}else{
			
			if(datalistfound[i]){
				# We need to transpose the data format, but to do that we first need to get it into a list:
				# If the data is like (var=val, ...) then we can just assign the value from eval(parse())
				# If the data is like (var <- val, ...) then the name of var is lost but var=val will be visible in the environment
				# So we need to use the environment to provide names for the variables with missing names				
				
				if(!grepl('list(', maindata[i], fixed=TRUE)){   # Need to make sure there is a list otherwise the templist won't be a list
					maindata[i] <- paste('list(\n', maindata[i], '\n)\n', sep='')
				}
				maindata[i] <- gsub(';',',',maindata[i],fixed=TRUE)   # semicolons will break the list	
				s <- try({
					tempenv <- new.env()
					templist <- eval(parse(text=maindata[i]), envir=tempenv)
					if(is.null(names(templist))){
						# If there are no names in g, all of the variables were assigned using <- so must be in tempenv:
						stopifnot(length(templist)==length(as.list(tempenv)))
						templist <- as.list(tempenv, all.names=TRUE)
					}
					# If some names are missing, only these will be in tempenv - so we could find the matches but it is better to just stop
					if(any(names(templist)=='')){
						stop('Some variable names were missing after attempting to parse the BUGS format list', call.=FALSE)
					}
				}, silent=TRUE)
				if(class(s)=='try-error')
					stop(paste('The following error was obtained while attempting to parse the list of ', if(data.type) 'data' else 'inits',  ':\n', as.character(s), sep=''), call.=FALSE)				
			}else{
				# We will use parse directly to get each of these and then turn them back into chars:
				s <- try({
					tempenv <- new.env()
					eval(parse(text=maindata[i]), envir=tempenv)
					templist <- as.list(tempenv, all.names=TRUE)
				}, silent=TRUE)
				if(class(s)=='try-error')
					stop(paste('The following error was obtained while attempting to parse the ', if(data.type) 'data' else 'inits', ':\n', as.character(s), sep=''), call.=FALSE)
				
			}
			
			# bugsdata is not necessarily the same as datalistfound[i] - it can be overridden
			if(bugsdata){
				#if(runjags.getOption('debug'))
				if(!transwarngiven)
					swcat('Note:  Transposing BUGS ', if(data.type) 'data' else 'inits', ' into R format\n')
				transwarngiven <- TRUE
				
				s <- try({				
	            	temp <- lapply(templist, function(x){
					if(is.null(dim(x)))
						return(x)

					dim(x) <- rev(dim(x))
					x <- aperm(x)
					return(x)
					})
				})
				if(class(s)=='try-error')
					stop(paste('The following unexpected error occured while converting the BUGS format ', if(data.type) 'data' else 'inits', ' into R format:\n', as.character(s), sep=''), call.=FALSE)

				templist <- temp
			}
			# If it is JAGS format data we can leave it alone		
			
			s <- try(maindata[i] <- dump.format(templist, checkvalid=FALSE), silent=TRUE)
			if(class(s)=='try-error')
				stop(paste('The following unexpected error occured while converting the ', if(data.type) 'data' else 'inits', ' back into dump format:\n', as.character(s), sep=''), call.=FALSE)
		}
		
	}
	
	
	# Get rid of any data that is empty or just whitespace:
	maindata <- maindata[maindata!='' & !grepl('^[[:space:]]*$', maindata)]
	if(length(maindata)==0)
		maindata <- ''
	
	if(length(maindata) > 1 && data.type){
		if(runjags.getOption('blockcombine.warning')) warning("More than 1 data block was found in the file.  Blocks were combined.", call.=FALSE)
		maindata <- paste(maindata, collapse="")
	}
	
	return(list(fixed=maindata, model=modelret))
}

#### Keeping in case of problems
oldsortjagsvsbugs <- function(maindata, bugs.format=bugs.format, data.type=TRUE){
	
	# This is all worded for data, but it applies equally to inits as well
	feedtxt <- if(data.type) 'A data' else 'An inits'
	
	# Note that the list( ) will already have been removed - datalistfound tells us if it was there:	
	datalistfound <- maindata[[2]]
	maindata <- maindata[[1]]
	
	modelret <- ''
	for(i in 1:length(maindata)){
		
		if(maindata[i]=='')
			next
		
		# Manual override mechanism for data types
    	tofinds <- tolower(maindata[i])
		bugsdata <- grepl("#bugsdata#", tofinds)   	# =fixed
		jagsdata <- grepl("#Rdata#", tofinds)  	# =fixed
		modeldata <- grepl("#modeldata#", tofinds)
		
    
		tdata <- maindata[i]
		
		browser()
		
		###########   DO SOME CONVERSIONS
		
		# Removes any  character between '#' and a new line
		# For regexp matching - '.' means any character, * means match the preceeding character any number of times

		# This would be greedy:
		#tdata <- gsub("#.*\n", "\n", maindata[i])
		# This is non-greedy matching:
		tdata <- gsub("#(.*?)\n", "\n", tdata, perl=TRUE)

		# This code will remove any "\n" from between "(" and ")" so that we can spot inappropriate commas later:
		# [^\\)] means any character except ) - it's double escaped for some reason
		greps <- gregexpr("\\([^\\)]*\n[^\\)]*\\)", tdata)
		s <- greps[[1]]
		e <- (s-1)+ attr(greps[[1]], "match.length")
		if(s[1]!=-1){ # If it doesn't find anything, s=-1
	
			finalstring <- strsplit(tdata, "", fixed=TRUE)[[1]]
			for(j in 1:length(s)){
	
				tstring <- finalstring[s[j]:e[j]]
				newstring <- gsub("\n", tstring)
				finalstring[s[j]:e[j]] <- newstring
	
			}
			tdata <- paste(finalstring, collapse="")
		}

		# Replaces any '),' with ')\n' (allowing for spaces in between) so that multiple vars on 1 line are allowed:
		tdata <- gsub(") *?, ", ")\n", tdata)

		
#
#		# Was doing helpful <- to = conversions but ditch them (you can't have .Dim <- structure or variable = value but that's up to the user)
#		temp <- strsplit(tdata[i], "")[[1]]
#		numbers <- which(temp=="=")
#		if(length(numbers)>0){
#			for(k in 1:length(numbers)){
#				tstring <- character(length=10)
#				for(j in 1:10){
#					tstring[j] <- paste(temp[pmax((numbers[k]-3-j):(numbers[k]-j), 1)], collapse="")
#				}	
#				if(all(tstring!=".Dim")) temp[numbers[k]] <- "<-"
#			}
#		}
#		tdata[i] <- paste(temp, collapse="")
#
		
		# Remove extra spaces from previously indented things:
		tdata <- gsub(",[[:space:]]*", ", ", tdata)
		
		# This code will change any "," to "\n" from between "<-" and "<-" unless a "(" or ")" are present:
		# [^\\)^\\(] means any character except ) or (
		# *? and perl makes it non-greedy matching
		# Need a while loop as the same character can't be used as the end point of one search and the start of another, so variable=value, variabile=value, variable=value would only remove the first (although variable=value, variable=value\nvariable=value, variable=value would be fine)
		s <- 0
		while(s!=-1){
			greps <- gregexpr("<-([^\\)^\\(]*?),([^\\)^\\(]*?)<-", tdata, perl=TRUE)
			s <- greps[[1]]
			e <- (s-1)+ attr(greps[[1]], "match.length")
			if(s[1]!=-1){ # If it doesn't find anything, s=-1
				finalstring <- strsplit(tdata, "", fixed=TRUE)[[1]]
				for(j in 1:length(s)){
	
					tstring <- finalstring[s[j]:e[j]]
					newstring <- gsub(",", "\n", tstring)
					finalstring[s[j]:e[j]] <- newstring
	
				}
				tdata <- paste(finalstring, collapse="")
			}
		}
		### Then make sure any new leading white space is removed (gsub is vectorised):
		tdata <- gsub("\n[[:space:]]*", "\n", tdata)
		####

		# As a final step, replace the removed , .Dim from data:
#		tdata <- gsub('\n.Dim', ',\n.Dim', tdata, fixed=TRUE)
		# And remove leading and trailing white space:
		tdata <- gsub('^ *\n', '', tdata)
		tdata <- gsub('\n *\n$', '\n', tdata)
		
		maindata[i] <- tdata
		###########   END OF CONVERSIONS
		
		# Fixed means it is to be read as data (jags or bugs column/row major order), model means part of the model:
		# models don't contain .Dim, c(, ' or " (I think)
		likelyfixed <- any(c(grepl(".Dim", maindata[i]), datalistfound[i], grepl("c\\(", maindata[i]), grepl("'", maindata[i]), grepl('"', maindata[i])))
		# fixed data doesn't contain [ or for(
		likelymodel <- data.type && any(c(grepl("\\[", maindata[i]), grepl("for *\\(", maindata[i])))
		# Since I'm removing comments I think it's safe to assume anything with = is fixeddata?, anything with .Dim is fixeddata, anything with list( is fixeddata, anything with c( is fixeddata, anything with [ is modeldata, anything with for( is modeldata.  Check all these and if we get a disparity then a warning?
    # Manual overrides are done above		

		if(modeldata && !data.type){
			warning('Ignoring an invalid #modeldata# specification in the initial values - assuming #Rdata#', call.=FALSE)
			modeldata <- FALSE
			bugsdata <- FALSE
			jagsdata <- TRUE
		}
		
		# If we don't have an override:
		if(!any(bugsdata, jagsdata, modeldata)){
			if((likelyfixed && likelymodel) || (!likelyfixed && !likelymodel)){
				if(runjags.getOption('debug'))
					swcat(feedtxt, ' type was not recognised - assuming fixed type\n', sep='')

				#print("try parsing here?  Have an option to default data type if unkown?")
				#for(i in 1:length(maindata)){
				#	err <- try(eval(parse(text=maindata[i])))
				#}

				# For now assume that default is jags type data:
				modeldata <- FALSE
				bugsdata <- bugs.format
				jagsdata <- !bugs.format
			}else{
				modeldata <- likelymodel
				bugsdata <- !likelymodel && bugs.format
				jagsdata <- !likelymodel && !bugs.format
			}
		}

		if(datalistfound[i] && !bugsdata)
			warning(paste(feedtxt, ' block that looks like a WinBUGS format was found - did you mean to specify bugs.format=TRUE?', sep=''), call.=FALSE)

		if(modeldata){
			if(modelret!='')
				stop('Multiple model style data blocks were found in the model - see the help file for ?read.jagsfile for options for specifying data', call.=FALSE)
			# Paste this data block onto the model and remove it from the data:
			modelret <- paste('data{\n', maindata[i], '\n}\n\n', sep='')
			maindata[i] <- ''
		}else{
			if(bugsdata){
				if(runjags.getOption('debug'))
					swcat('Transposing BUGS data\n')
				
				# We need to transpose the data format:
				s <- try({
				temp <- list.format(maindata[i], checkvalid=FALSE)
	            temp <- lapply(temp, function(x){
					if(is.null(dim(x)))
						return(x)

					dim(x) <- rev(dim(x))
					x <- aperm(x)
					return(x)
					})
				maindata[i] <- dump.format(temp, checkvalid=FALSE)
				})
				if(class(s)=='try-error')
					stop('An unexpected error occured while converting the BUGS format data into JAGS format', call.=FALSE)					
			}
			# If it is JAGS format data we can leave it alone				
		}
		
	}
	
	
	# Get rid of any data that is empty or just whitespace:
	maindata <- maindata[maindata!='' & !grepl('^[[:space:]]*$', maindata)]
	if(length(maindata)==0)
		maindata <- ''
	
	if(length(maindata) > 1 && data.type){
		if(runjags.getOption('blockcombine.warning')) warning("More than 1 data block was found in the file.  Blocks were combined.", call.=FALSE)
		maindata <- paste(maindata, collapse="")
	}
	
	return(list(fixed=maindata, model=modelret))
}


winbugs.extract.small <- function(find, string){
  
  if(identical(string, '\n') || identical(string, '\r') || identical(string, ''))
    return('')
  
  # Eliminate \t and add the # - allow 0+ spaces:
  string <- gsub('\t',' ',string)
  grfind <- paste('# *',find,' *#',sep='')
	# Separate into lines:
  split <- strsplit(string,split='[\n\r]')[[1]]
  # Eliminate lines without the target:
  split <- split[grepl(grfind,split,ignore.case=TRUE)]
  
  # If none, return:
  if(length(split)==0)
    return('')
  
  # Now remove anything before the target identifier:
  split <- gsub(paste('^.*',grfind,sep=''),'',split,ignore.case=TRUE)
  
  # Now just separate on the allowed characters and delist:
  found <- unlist(strsplit(split, '[\\,,;,&]'))
  
  # Remove trailing and leading spaces:
  found <- gsub('^ *', '', found)
  found <- gsub(' *$', '', found)
  
  # Stop if an illegal space is found unless this is modules or factories:
  if(any(grepl(' ', found)) && !find%in%c('modules','factories')){
    problem <- which(grepl(' ', found))
    stop(paste('Invalid space in the following specification(s) for #', find, '#: ', paste(found[problem], collapse=', ')), call.=FALSE)
  }
  
  # And return:
  return(found)
  
}


# Keep this in case it's needed - old and slow (but more general?) version
old.winbugs.extract.small <- function(find, string){

split <- strsplit(string, "")[[1]]

newstring <- ""

newlinelast = found <- FALSE
find.no = hash <- 0

for(i in 1:length(split)){
		
	if(found){
		
		if(split[i]=="#"){
			hash <- hash + 1
			next
		}

		if(hash!=2) next
		
		if(any(split[i]==c(",", ";", ":", "&"))){
			
			find.no <- find.no + 1
			next
		}
		
		if(is.na(newstring[find.no])) newstring[find.no] <- ""
		
		if(any(split[i]==c("", " ", "\t", "@", "%"))) next
		if(any(split[i]==c("\n", "\r"))){
			found <- FALSE
		}else{
			newstring[find.no] <- paste(newstring[find.no], split[i], sep="")
		}
				
	}

	if(!all(split[i]==c("", " ", "\t"))){
		temp <- split[i:length(split)]
		temp <- temp[temp!=" " & temp!= "" & temp!="\t"]
		if(paste(temp[1:(length(strsplit(find, "")[[1]]))], collapse="") == find){ # newlinelast not necessary for extract.small
			found <- TRUE
			hash <- 1
			split[min(which(split=="#")[which(split=="#")>=i])] <- ""
			find.no <- find.no + 1
		}
	
	
		if(any(split[i]==c("\n", "\r"))){
			newlinelast <- TRUE
		}else{
			newlinelast <- FALSE
			}
	}

}
for(i in 1:length(newstring)){
	temp <- strsplit(newstring[i], "")[[1]]
	temp[temp=="="] <- "<-"
	newstring[i] <- paste(temp, collapse="")
}
return(newstring[newstring!="" & newstring!=" "])
}


normalise.mcmcfun <- function(mcmc.list, normalise = TRUE, warn = TRUE, remove.nonstochastic = TRUE){

	if(class(mcmc.list)=="mcmc") mcmc <- mcmc.list(mcmc.list) else mcmc <- mcmc.list

	if(class(mcmc)!="mcmc.list") stop("Object to be normalised must be an mcmc list or mcmc object")
		
	vnames <- lapply(mcmc, dimnames)
	
	# Special case of 0 iterations:
	if(niter(mcmc)==0){
		truestochastic <- rep(TRUE, nvar(mcmc))
		names(truestochastic) <- vnames[[1]][[2]]
		nonstochastic=semistochastic <- !truestochastic
		
	}else{
	
		# ALWAYS cehck stochastic, but sometimes don't remove them
		variances <- sapply(mcmc, function(x){			
			return(apply(x,2,function(y){
			
				usey <- y[y!=Inf & y!=-Inf & !is.na(y)]
				if(length(usey)<=1) v <- 0 else v <- var(usey)
				return(v)
			
			}))}, simplify='array')
	
		dim(variances) <- c(nvar(mcmc),nchain(mcmc))
		dimnames(variances) <- list(vnames[[1]][[2]], NULL)
		anyvariancezero <- variances==0
	
		means <- sapply(mcmc, function(x){			
			return(apply(x,2,function(y){
			
				usey <- y[y!=Inf & y!=-Inf & !is.na(y)]
				if(length(usey)==0) v <- NA else v <- mean(usey)
				return(v)
			
			}))}, simplify='array')
	
		dim(means) <- c(nvar(mcmc),nchain(mcmc))
		dimnames(means) <- list(vnames[[1]][[2]], NULL)
	
		# If there is only one chain then meansdiffer is just taken from variance within the chain i.e. ignored:
		meansdiffer <- apply(means,1,var)>0
		if(nchain(mcmc)==1)
			meansdiffer[] <- !apply(anyvariancezero,1,any)
	
		# There are 3 options:
		# variances > 0 between chains but == 0 within chains ->  semi-stochastic
		# variance == 0 in some but not all chains -> semi-stochastic
		# variances == 0 in all chains but means differ between -> semi-stochastic
		# variances == 0 in all chains and means identical -> non-stochastic
		# all variances > 0 -> stochastic
		# 1 iteration variance will be 0 so counted as semi-stochastic
	
		# 1 iteration and chain -> everything non-stochastic
		if(nchain(mcmc)==1 && niter(mcmc)==1){
			anyvariancezero[] <- TRUE
			meansdiffer[] <- FALSE
		}
	
		nonstochastic <- apply(anyvariancezero,1,all) & !meansdiffer
		removed <- names(nonstochastic)[nonstochastic]
		if(length(removed)>0 && niter(mcmc)>1){
			warnmessage <- paste("Note: The monitored variable", if(sum(nonstochastic)>1) "s", " '", if(sum(nonstochastic)>1) paste(removed[1:(length(removed)-1)], collapse="', '"), if(sum(nonstochastic)>1) "' and '", removed[length(removed)], "' appear", if(sum(nonstochastic)==1) "s", " to be non-stochastic; ", if(sum(nonstochastic)>1) "they" else "it", " will not be included in the convergence diagnostic", sep="")	
			if(warn==TRUE)
				swcat(warnmessage,"\n")
			if(warn=="warning")
				warning(warnmessage, call.=FALSE)
		}
	
		truestochastic <- apply(!anyvariancezero,1,all)  # NOT necessary for the means to be different - they might be the same by chance
		
		#incomplete?:
		#semistochastic <- (apply(anyvariancezero,1,any) != apply(anyvariancezero,1,all)) | (apply(anyvariancezero,1,all) & meansdiffer)
		semistochastic <- !(nonstochastic | truestochastic)
		
		removed <- sort(names(semistochastic)[semistochastic])
		if(length(removed)>0 && niter(mcmc)>1){
			warnmessage <- paste("Note: The monitored variable", if(sum(semistochastic)>1) "s", " '", if(sum(semistochastic)>1) paste(removed[1:(length(removed)-1)], collapse="', '"), if(sum(semistochastic)>1) "' and '", removed[length(removed)], "' appear", if(sum(semistochastic)==1) "s", " to be stochastic in one chain but non-stochastic in another chain; ", if(sum(semistochastic)>1) "they" else "it", " will not be included in the convergence diagnostic", sep="")	
			if(warn==TRUE)
				swcat(warnmessage,"\n")
			if(warn=="warning")
				warning(warnmessage, call.=FALSE)
		}
	
		if(sum(semistochastic)+sum(truestochastic)+sum(nonstochastic) != length(vnames[[1]][[2]])){
			swcat('True stochastic:\n')
			print(truestochastic)
			swcat('Semi stochastic:\n')
			print(semistochastic)
			swcat('Non stochastic:\n')
			print(nonstochastic)
			swcat('Variable names:\n')
			print(vnames[[1]][[2]])
			if(runjags.getOption('debug'))
				browser()
			
			stop('An unexpected error occured determining which variables were and were not stochasic - please file an error report to the package maintainer', call.=FALSE)			
		}
	}
	
#		if(all(nontruestochastic))
#			stop("All monitored variables appear to be non-stochastic; try adding monitored variables either using #monitor# in the model code or by specifying a monitor argument to run.jags", call.=FALSE)

	# remove.nonstochastic is STRICT nonstochastic only (i.e. semistochastic and stochastic both left alone)
	if(remove.nonstochastic && any(nonstochastic))
		mcmc <- mcmc[,-which(nonstochastic),drop=FALSE]
		
	
	if(normalise && !all(nonstochastic)){
		if(niter(mcmc)>1000) use <- sample(1:niter(mcmc), size=1000, replace=FALSE) else use <- 1:niter(mcmc)

		shap.res <- apply(combine.mcmc(mcmc, collapse.chains=TRUE, add.mutate=FALSE), 2, function(x){
			rv <- 1
			if(!any(x==Inf | x==-Inf | is.na(x)) && var(x)!=0){
				# We might get an error if we sample 1000 funny values - in which case leave them alone (likely to be small variance so transform probably unnecessary anyway)
				suppressWarnings(success <- try({
					if(all(x > 0)){
						if(all(x < 1)){
							if(stats::shapiro.test(x[use])$p.value < stats::shapiro.test(log(x[use]/(1-x[use])))$p.value) rv <- 3
						}else{
							if(stats::shapiro.test(x[use])$p.value < stats::shapiro.test(log(x[use]))$p.value) rv <- 2
						}
					}
				}, silent=TRUE))
			}
			return(rv)
		})

		change <- which(shap.res!=1)
		for(parameter in change){
			for(chain in 1:nchain(mcmc)){
				newvalues <- unlist(mcmc[[chain]][,parameter])
				if(shap.res[parameter]==3) newvalues <- log(newvalues/(1-newvalues)) else newvalues <- log(newvalues)
				mcmc[[chain]][,parameter] <- newvalues
			}
		}
	}
	for(i in 1:length(mcmc))
		dimnames(mcmc[[i]])[[1]] <- vnames[[i]][[1]]

	if(class(mcmc.list)=="mcmc"){
		return(list(mcmc=mcmc[[1]], truestochastic=truestochastic, semistochastic=semistochastic, nonstochastic=nonstochastic))
	}else{
		return(list(mcmc=mcmc, truestochastic=truestochastic, semistochastic=semistochastic, nonstochastic=nonstochastic))
	}
}

safe.autocorr.diag <- function(x, ...){
	if(niter(x)==1){
		y <- matrix(NA, nrow=1, ncol=nvar(x), dimnames=list('Lag 1', dimnames(x[[1]])[[2]]))
		return(y)
	}
		
	y <- autocorr.diag(x[,1],...)
	if(nvar(x)>1) for(i in 2:nvar(x)) y <- cbind(y, autocorr.diag(x[,i],...))
	dimnames(y)[[2]] <- dimnames(x[[1]])[[2]]
	return(y)
}

safe.gelman.diag <- function(x, warn=TRUE,...){
	
	success <- try(gelman <- gelman.diag(x, ...), silent=TRUE)
	if(inherits(success, 'try-error')){
		
		nvars <- nvar(x)
		psrfs <- matrix(ncol=2, nrow=nvars, dimnames=list(varnames(x), c("Point est.", "97.5% quantile")))
		
		success <- try({
			for(i in 1:nvars){
			psrfs[i,] <- gelman.diag(x[,i,drop=FALSE], ...)$psrf
		}
		}, silent=TRUE)
		if(inherits(success, 'try-error')){
			#name <- new_unique("gelman.failed", ".Rsave")
			#save(failedmcmc, file=name)
			stop("An error occured while calculating the Gelman-Rubin statistic")
		}
		
		if(warn) swcat("Note:  Unable to calculate the multivariate psrf\n")
		
		y <- list(psrf=psrfs, mpsrf="Unable to calculate multivariate psrf")
		
		class(y) <- "gelman.diag"
		return(y)
	}else{
		return(gelman)
	}
	
}

getargs <- function(functions, passed, returnall=TRUE, otherfnames=character(0)){
	
	N <- length(functions)
	args <- vector('list', length=N)
	names <- vector('list', length=N)
		
	# Argument names for the functions relative to this function:
	for(i in 1:N){
		args[[i]] <- as.list(formals(get(functions[i], sys.frame(sys.parent(n=1)))))
		args[[i]] <- args[[i]][names(args[[i]])!="..."]
		names[[i]] <- names(args[[i]])
	}
	
 	argnames <- unique(unlist(names))
 	argmatch <- pmatch(names(passed), argnames)

 	if(any(is.na(argmatch))){
     
     nomatches <- names(passed)[which(is.na(argmatch))]
 		
		functions <- c(otherfnames, functions)
 		functstring <- paste(if(length(functions)>1) paste(functions[1:(length(functions)-1)], collapse="', '"), if(length(functions)>1)"' or '", functions[length(functions)], "' function", if(length(functions)>1) "s", sep="")
 		
 		argstring <- paste(if(length(nomatches)>1) "s", " '", if(length(nomatches)>1) paste(nomatches[1:(length(nomatches)-1)], collapse="', '"), if(length(nomatches)>1)"' or '", nomatches[length(nomatches)], "'", sep="")

		stop(paste("unused argument(s)", argstring, " (no unambiguous match in the '",functstring,")", sep=""), call.=FALSE)
 	}
 	
	names(passed) <- argnames[argmatch]
	passed <- passed[!is.na(argmatch)]
 	
	if(returnall){
		# Now get defaults from specified functions, giving priority to earlier functions if arguments appear in more than 1:
		alreadymatched <- names(passed)
		for(i in 1:N){
			newget <- names[[i]][!(names[[i]] %in% alreadymatched)]
			newargs <- lapply(newget, function(x) try(as.expression(get(x, pos=args[[i]])), silent=TRUE))
			names(newargs) <- newget
			newargs <- newargs[!sapply(newargs,class)=="try-error"]
			passed <- c(passed, newargs)
			alreadymatched <- c(alreadymatched, newget)
		}
	}
	return(passed)
}

tailf <- function(file, start=1, refresh=0.1, min.static=1, max.static=Inf, stop.function=function() return(FALSE), stop.text=character(0), print=TRUE, return=!print){

	done <- FALSE

	readupto <- 1

	# Allow the simulation to start:
	if(!file.exists(file)) Sys.sleep(1)

	tryCatch({
	if(!file.exists(file)) stop("The named file does not exist")

	linesto <- start+1
	going <- TRUE
	static <- 0
	lastline <- ""

	text <- ""

	# Wait for file to start being written to:
	repeat{
		suppressWarnings(output <- readLines(file))
		# Catch occasional error with na being introduced:
		if(is.na(lastline)) lastline <- ""
		if(length(output)==0){
			static <- static+1
			if(static > max.static) break
		}else{
			break
		}
	}

	repeat{
	Sys.sleep(refresh)

	suppressWarnings(output <- readLines(file))
	# Catch occasional error with na being introduced:
	if(is.na(lastline)) lastline <- ""

	if(length(output)==(linesto-1)){
		if(output[length(output)]==lastline){
			static <- static+1
			if(static > max.static) break
		}else{
			static <- 0
			newsameline <- output[readupto]		
			new <- strsplit(newsameline, '')[[1]]
			last <- length(strsplit(lastline, '')[[1]])
			new <- paste(new[(last+1):length(new)], collapse='')
			if(print) cat(new)
			text <- paste(text, new, sep='')
			lastline <- output[length(output)]
		}
	}else{
		if(output[linesto-1]!=lastline){
			new <- output[linesto-1]
			new <- strsplit(new, '')[[1]]
			last <- length(strsplit(lastline, '')[[1]])
			new <- paste(new[(last+1):length(new)], collapse='')
			if(print) cat(new)
			text <- paste(text, new, sep='')
		}
		static <- 0
		new <- paste(output[linesto:length(output)], collapse='\n')
		if(print) cat('\n', new, sep='')
		text <- paste(text, '\n', new, sep='')
		linesto <- length(output)+1
		lastline <- output[length(output)]
	}

	flush.console()
	if(static > min.static & stop.function()) break
	
	allnewoutput <- paste(output[readupto:length(output)],collapse="\n")
	readupto <- length(output)
	foundtextmatch <- any(sapply(as.list(stop.text), grepl, x=allnewoutput))
	if(foundtextmatch) break

	}

	done <- TRUE

	}, finally={
	
		retval <- list(text=text, lines=length(output), interrupt=!done)
		if(return) return(retval)
	
		})
}

prettifytable <- function(x, digits=5, colsequal=FALSE, nastring="", psrfcoldollar=FALSE){

	formatted <- formatC(x, format="fg", digits=digits, width=-1)
	absx <- abs(x)
	absx[is.na(absx)] <- 1
	formatted[absx<10^-digits] <- formatC(x[absx<10^-digits], digits=digits, width=-1)
	
	formatted <- gsub("NA",nastring,formatted)
	formatted <- gsub("NaN",nastring,formatted)
	if(any(psrfcoldollar)){
		stopifnot(length(psrfcoldollar)==nrow(formatted))		
#		formatted[psrfcoldollar,'psrf'] <- '$'
		formatted[psrfcoldollar,'psrf'] <- paste(formatted[psrfcoldollar,'psrf'], ' $', sep='')
	}
	
	# Put column names on as well:
	formatted <- rbind(dimnames(formatted)[[2]],formatted)
	dimnames(formatted) <- list(dimnames(formatted)[[1]],rep("",ncol(formatted)))
	
	if(colsequal){
		retval <- format(formatted, justify="right")
	}else{
		retval <- apply(formatted, 2, format, justify="right")
	}	
	
	# noquote suppresses the "" from the character string:
	return(noquote(retval))
	
}

checkmodfact <- function(tocheck, type){
	if(identical(tocheck,'')) tocheck <- list()
	stopifnot(type%in%c('module','factory'))
  	
  	# In case any are blank:
  	if(length(tocheck)>0)
		tocheck <- tocheck[sapply(tocheck, function(x) return(!all(x=='')))]
  	  
	nl <- switch(type, module=2, factory=3)

  	if(is.character(tocheck)){
		if(length(tocheck)>0 && any(grepl(',', tocheck, fixed=TRUE)))
			stop('Use of commas in module or factory specifications is not allowed - separate name type and status with a space', call.=FALSE)
		tocheck <- gsub('[\\(\\)]',' ',tocheck)
	  	tocheck <- strsplit(gsub('[[:space:]]+', ' ', tocheck),' ')
  	}
  	
	if(!is.list(tocheck)){
		stop(paste('Invalid ', type, ' specification - it must be either a character vector or a list', sep=''), call.=FALSE)
	}
	
	# If a blank list return '':
	if(identical(list(), tocheck)) return('')
	
  	validated <- lapply(tocheck, function(x){
		origx <- x

	  	# If not specified, assume it wants to be on:
	  	if(length(x) < nl)
			x <- c(x, "TRUE")
		
		# Check length is correct:
		if(length(x)!=nl)
			stop(paste('Incorrect number of elements for ', type, ' specification "', paste(origx,collapse=' '), '": ', length(origx), ' found but ', nl, ' expected', sep=''), call.=FALSE)
		
	  	# Replace on with TRUE and off with FALSE:
	  	x[nl] <- gsub('on','TRUE',x[nl])
	  	x[nl] <- gsub('ON','TRUE',x[nl])
	  	x[nl] <- gsub('off','FALSE',x[nl])
	  	x[nl] <- gsub('OFF','FALSE',x[nl])
	  	x[nl] <- as.character(as.logical(x[nl]))
		
		# Check factory type is valid:
	  	if(type=='factory' && !x[2]%in%c('sampler','monitor','rng'))
			stop(paste('The type of ', type, ' specification "', paste(origx,collapse=' '), '" must be one of sampler, monitor or rng', sep=''), call.=FALSE)
	
	  	# Check they are all logicable:
	  	if(is.na(x[nl]))
	  		stop(paste('The status "', origx[nl], '" of ', type, ' specification "', paste(origx, collapse=' '), '" is not interpretable as logical', sep=''), call.=FALSE)
	
		return(x)
	})

  validated <- unique(validated)
  if(type=='module')
    tocheck <- sapply(validated,function(x) return(x[1]))
  if(type=='factory')
    tocheck <- sapply(validated,function(x) return(paste(x[1],x[2],sep=' ')))
	if(length(unique(tocheck))!=length(tocheck))
		stop(paste('Replicated ', type, ' name(s) with conflicting status in "', paste(sapply(validated,paste,collapse=' '), collapse=', '), '"', sep=''), call.=FALSE)	
	if(identical(list(), validated)) validated <- ''
	
	return(validated)
}

checkvalidrunjagsobject <- function(runjags.object){

	if(!is.runjags(runjags.object))
		stop("The output of a runjags function must be supplied", call.=FALSE)
	
	# Add some info that may not be available from old (saved) runjags objects:

	# Only add things if not debugging:
	if(!runjags.getOption('debug')){
		
	  # Added for version 2.0:
		if(is.null(runjags.object$summary.pars))
			runjags.object$summary.pars <- getdefaultsummarypars()
		if(is.null(runjags.object$summary.available))
			runjags.object$summary.available <- FALSE
		if(is.null(runjags.object$deviance.table))
			runjags.object$deviance.table <- NA
		if(is.null(runjags.object$deviance.sum))
			runjags.object$deviance.sum <- NA
		if(is.null(runjags.object$pd))
			runjags.object$pd <- NA
		if(is.null(runjags.object$samplers))
			runjags.object$samplers <- NA
		if(is.null(runjags.object$semistochastic))
			runjags.object$semistochastic <- rep(nvar(runjags.object$mcmc), FALSE)
		if(is.null(runjags.object$truestochastic))
			runjags.object$truestochastic <- runjags.object$semistochastic
	  	runjags.object$nonstochastic=!(runjags.object$truestochastic | runjags.object$semistochastic)
	  	if(is.list(runjags.object$dic) && is.null(runjags.object$dic$meanpopt))
			runjags.object$dic$meanpopt <- NA
		if(is.null(runjags.object$response))
			runjags.object$response <- NA
		if(is.null(runjags.object$residual))
			runjags.object$residual <- NA
		if(is.null(runjags.object$fitted))
			runjags.object$fitted <- NA
	
		# Not the best way of doing it - add.summary doesn't update the version:
	#	if(numeric_version(runjags.object$runjags.version[1]) < 2){
	#		runjags.object$summaries <- 'Summary statistics not available - see ?add.summary'
	#	}
		
		defmo <- getdefaultmethodoptions()
		if(is.null(runjags.object$method.options)){
			runjags.object$method.options <- defmo
		}

		if(is.null(runjags.object$method.options$n.sims))
			runjags.object$method.options$n.sims <- defmo$n.sims
		if(is.null(runjags.object$method.options$cl))
			runjags.object$method.options$cl <- defmo$cl
		if(is.null(runjags.object$method.options$remote.jags))
			runjags.object$method.options$remote.jags <- defmo$remote.jags
		if(is.null(runjags.object$method.options$by))
			runjags.object$method.options$by <- defmo$by
	  	if(is.null(runjags.object$method.options$progress.bar))
	    	runjags.object$method.options$progress.bar <- defmo$progress.bar

		if(!all(c('n.sims', 'cl', 'remote.jags', 'by', 'progress.bar', 'jags', 'silent.jags', 'jags.refresh', 'batch.jags')%in%names(runjags.object$method.options))){
			runjags.object$method.options <- defmo
		}

	  
	}else{
		if(!all(c('plots','vars','mutate','psrf.target','normalise.mcmc','modeest.opts','confidence','autocorr.lags','custom','silent.jags','plots','plot.type','col','summary.iters','trace.iters','separate.chains','trace.options','density.options','histogram.options','ecdfplot.options','acplot.options')%in%names(runjags.object$summary.pars) || length(runjags.object$summary.pars)!=21)){
			cat('Invalid summary.pars\n')
			browser()
		}
	}
	
	# Less strict - don't check for the presence of summary statistics as they may not have been added yet:
	if(!all(c('mcmc', 'end.state', 'burnin', 'sample', 'thin', 'model', 'data', 'monitor', 'noread.monitor', 'modules', 'factories', 'response', 'residual', 'fitted', 'method', 'method.options', 'timetaken', 'summary.pars', 'deviance.table', 'deviance.sum') %in% names(runjags.object)))
		stop("Invalid runjags.object provided; the output of a runjags function (with class 'runjags') must be supplied", call.=FALSE)
	
  # rjags is NOT required as it is sometimes removed to decompile it
	if(!all(c('n.sims', 'cl', 'remote.jags', 'by', 'progress.bar', 'jags', 'silent.jags', 'jags.refresh', 'batch.jags')%in%names(runjags.object$method.options)))
		stop("Invalid runjags.object provided (invalid method options); the output of a runjags function (with class 'runjags') must be supplied", call.=FALSE)
	
	
	if(length(runjags.object$response)!=1 || (!is.na(runjags.object$response) && !is.character(runjags.object$response)))
		stop('The response variable (if supplied) must be a single character variable')
	if(length(runjags.object$residual)!=1 || (!is.na(runjags.object$residual) && !is.character(runjags.object$residual)))
		stop('The residual variable (if supplied) must be a single character variable')
	if(length(runjags.object$fitted)!=1 || (!is.na(runjags.object$fitted) && !is.character(runjags.object$fitted)))
		stop('The fitted variable (if supplied) must be a single character variable')
	
	runjags.object$factories <- checkmodfact(runjags.object$factories, 'factory')
	runjags.object$modules <- checkmodfact(runjags.object$modules, 'module')
		
	invisible(runjags.object)
}

checkvalidmonitorname <- function(monitor){
	
	if(is.logical(monitor) && length(monitor)>0 && !identical(NA, monitor)){
		if(any(is.na(monitor))){
			warning('Setting NA values in specified logical argument for variable selection to FALSE', call.=FALSE)
			monitor[is.na(monitor)] <- FALSE
		}
		return(monitor)
	}
	
	monitor <- as.character(monitor)
	
	# Remove spaces from monitors
	monitor <- gsub('[[:space:]]','',monitor)
	
	# Don't allow comments in monitor names
	problem <- grepl('#', monitor)
	if(any(problem))
		stop(paste('Invalid monitor name(s) "', paste(monitor[problem],collapse='", "'), '" - comment (#) characters are not allowed', sep=''), call.=FALSE)
		
	# Look for ,, and [, and ,] and give an error
	problem <- grepl('[,', monitor, fixed=TRUE) | grepl(',]', monitor, fixed=TRUE) | grepl(',,', monitor, fixed=TRUE)
	if(any(problem))
		stop(paste('Invalid monitor name(s) "', paste(monitor[problem],collapse='", "'), '" - empty indexes are not allowed', sep=''), call.=FALSE)

	# Quotation marks are allowed (as are $ and ^), as we use them for matching specific variables:
#	problem <- grepl('"', monitor, fixed=TRUE) | grepl("'", monitor, fixed=TRUE)
#	if(any(problem))
#		stop(paste('Invalid monitor name(s): ', paste(monitor[problem],collapse=', '), ' - quotation marks are not allowed', sep=''), call.=FALSE)
	
	# Look for unmatched brackets:
	problem <- (grepl('[', monitor, fixed=TRUE) & !grepl("]", monitor, fixed=TRUE)) | (!grepl('[', monitor, fixed=TRUE) & grepl("]", monitor, fixed=TRUE))
	if(any(problem))
		stop(paste('Invalid monitor name(s): ', paste(monitor[problem],collapse=', '), ' - unmatched square bracket', sep=''), call.=FALSE)

	# Look for , without [:
	problem <- grepl(',', monitor, fixed=TRUE) & !grepl("[", monitor, fixed=TRUE)
	if(any(problem))
		stop(paste('Invalid monitor name(s): ', paste(monitor[problem],collapse=', '), ' - commas are only allowed within square brackets', sep=''), call.=FALSE)
	
	# The word 'to' is going to give a problem:
	if(!all(is.na(monitor)) && any(monitor=='to'))
		warning('Use of the monitor name "to" may cause problems with some JAGS methods', call.=FALSE)
	
	monitor <- expandindexnames(monitor)
	return(monitor)
}

# Required for modified read.coda function to work properly:
expandindexnames <- function(names){
	f <- function(x){
	  nums <- suppressWarnings(as.numeric(x))
	  keep <- x[is.na(nums)]
	  nums <- nums[!is.na(nums)]
	  npairs <- length(nums)/2
	  if(npairs>0){
	    indexes <- expand.grid(lapply(1:npairs,function(p){
	      s <- ((p-1)*2)+1
	      return(nums[s]:nums[s+1])
	    }))
	    keep <- t(matrix(keep, nrow=length(keep), ncol=nrow(indexes)))
	    for(i in 1:npairs) keep[,(i*2)] <- indexes[,i]
	    keep <- apply(keep,1,paste,collapse='')
	  }
	  return(keep)
	}
	
	# Mark indexes using : using hashes in a predetermined way:
	newnames <- gsub('([[:digit:]]+):([[:digit:]]+)', '##\\1#\\2#', names)			
	# Now split on comment char:
	split <- strsplit(newnames,'#')
	# And recombine with the indexes explicitly rolled out:
	toret <- as.character(unlist(lapply(split,f)))
	# Remove any spaces:
	toret <- gsub(' ','',toret)
	# Check the bracket isn't empty:
	if(any(grepl('[]',toret,fixed=TRUE))){
		stop(paste('Empty indexes provided for variable(s) ', paste(toret[grepl('[]',toret,fixed=TRUE)], collapse=', '), ' - ensure that the full range of indices are specified using a colon (e.g. var[1:2,1])', sep=''), call.=FALSE)
	}
	
	if(length(toret)==0 || identical(toret, '') || identical(toret, NA))
	  return(toret)
	
	if(grepl(',,',toret,fixed=TRUE) || grepl('[,',toret,fixed=TRUE) || grepl(',]',toret,fixed=TRUE))
		stop(paste('Ambiguous index entry provided for variable(s) ', paste(toret[grepl(',,',toret,fixed=TRUE) || grepl('[,',toret,fixed=TRUE) || grepl(',]',toret,fixed=TRUE)], collapse=', '), ' - ensure that all indices are specified using a colon (e.g. var[1:2,1])', sep=''), call.=FALSE)
	
	return(toret)	
}

matchvars <- function(vars, names, exactneeded=NA, testfound=TRUE){
	
	if(length(names)==0)
		stop('No valid variable names in the object supplied', call.=FALSE)
	
	if(is.logical(vars) && length(vars)>0){
		if(length(vars)!=length(names)){
			stop(paste("The length of the logical vector specified to 'vars' (", length(vars), ") does not match the number of monitored variables (", length(names), ")", sep=""), call.=FALSE)
		}
		if(any(is.na(vars)))
			stop("Missing values are not allowed in the logical vector specified to 'vars'", call.=FALSE)
		
		return(which(vars))
	}
	
	vars <- as.character(na.omit(vars))
	# vars is what to find, names is the MCMC object varnames to find it in
	
	if(length(vars)>0){
		
		vars <- expandindexnames(vars)
	#	matched <- vapply(vars, function(m) return(grepl(paste("^",m,sep=""),names)), logical(length(names)))
		
		matched <- vapply(vars, function(m) return(grepl(m,paste("^",names,"$",sep=""),fixed=TRUE)), logical(length(names)))	
		exact <- vapply(vars, function(m) return((gsub("'","",gsub('"','',m,fixed=TRUE),fixed=TRUE)) == names), logical(length(names)))
		namesnoindex <- sapply(strsplit(names,'[',fixed=TRUE), function(x) return(x[[1]]))
		nameswithoutindex <- vapply(vars, function(m) return((gsub("'","",gsub('"','',m,fixed=TRUE),fixed=TRUE)) == namesnoindex), logical(length(names)))
		varsnoindex <- sapply(strsplit(vars,'[',fixed=TRUE), function(x) return(x[[1]]))
		varswithoutindex <- vapply(varsnoindex, function(m) return((gsub("'","",gsub('"','',m,fixed=TRUE),fixed=TRUE)) == names), logical(length(names)))
		
		if(is.null(dim(matched))){
			dn <- names(matched)
			dim(matched) <- c(length(names), length(vars))
			dimnames(matched) <- list(NULL, dn)
		}
		if(is.null(dim(exact))){
			dn <- names(exact)
			dim(exact) <- c(length(names), length(vars))
			dimnames(exact) <- list(NULL, dn)
		}
		if(is.null(dim(nameswithoutindex))){
			dn <- names(nameswithoutindex)
			dim(nameswithoutindex) <- c(length(names), length(vars))
			dimnames(nameswithoutindex) <- list(NULL, dn)
		}
		if(is.null(dim(varswithoutindex))){
			dn <- names(varswithoutindex)
			dim(varswithoutindex) <- c(length(names), length(vars))
			dimnames(varswithoutindex) <- list(NULL, dn)
		}
		
		if(identical(exactneeded, NA))
			exactneeded <- t(matrix((grepl("'",vars,fixed=TRUE) | grepl('"',vars,fixed=TRUE)), ncol=length(names), nrow=length(vars)))
		
		if(!all(dim(matched)==dim(exact)) || !all(dim(matched)==dim(nameswithoutindex)) || !all(dim(matched)==dim(varswithoutindex))){
			print(matched)
			print(exact)
			stop('An unexpected error occured while matching variable names - please file a bug report to the package maintainer including the matrices printed above', call.=FALSE)
		}
		
		# This still needs to be vectorised like this in case exactneeded depends on var:
		combined <- (matched & !exactneeded) | (varswithoutindex & !exactneeded) | exact | nameswithoutindex
		selected <- lapply(1:ncol(combined), function(x){
			return(which(combined[,x]))
		})
		names(selected) <- dimnames(combined)[[2]]
		# selected <- apply( (matched & !exactneeded) | (exact & exactneeded) , 2, function(x) return(if(any(x)) which(x) else 0))

		if(any(sapply(selected,length)==0) && testfound)
			stop(paste("No matches found for the following variable name(s): ", paste(vars[sapply(selected,length)==0],collapse=','), sep=''), call.=FALSE)

		selected <- unique(unlist(selected))
		
	}else{
			
		selected <- 1:length(names)
	}
	
	return(selected)
}

checkvalidforjags <- function(object){
	
	if(length(object)==1 && is.na(object)) return(list(valid=TRUE, probstring=""))
	
	if(class(object)=="runjagsdata" || class(object)=="runjagsinits") class(object) <- "character"
	if(class(object)!="list" && class(object)!="character") return(list(valid=FALSE, probstring="object must be either a named list or a character vector in the R dump format"))
	
	if(!is.list(object)) object <- list.format(object, checkvalid=FALSE)
	
	if(any(names(object) == "")){
		return(list(valid=FALSE, probstring="missing variable name(s)"))
	}

	if(!length(unique(names(object))) == length(object)){
		return(list(valid=FALSE, probstring="duplicated variable name(s)"))
	}
	
	problems <- sapply(object, function(x){
		
		# Catch potential problems with the data being passed through:
		if(length(x)==0){
			return("")
		}
		if(is.null(x)){
			return("NULL")
		}
		if(inherits(x, "data.frame")){
			return("inherits from class 'data.frame' - try converting it to a valid type using as.matrix")
		}
		if(length(x)==0){
			return("length zero")
		}
		if(class(x)=="logical" && !all(is.na(x))){
			return("TRUE/FALSE")
		}
		if(class(x)=="character" && !all(is.na(x))){
			return("character")
		}
		if(class(x)=="factor" && !all(is.na(x))){
			return("factor")
		}
		if(any(x==Inf, na.rm=TRUE)){
			return("Inf")
		}
		if(any(x==-Inf, na.rm=TRUE)){
			return("Inf")
		}
		return("")
	})	
	
	problems[names(problems)==".RNG.name"] <- ""
	
	if(all(problems=="")){
		return(list(valid=TRUE, probstring=""))
	}else{
		problems <- problems[problems!=""]
		probstring <- paste("invalid variable value(s) - ", paste(names(problems), " (", problems, ")", sep=""), collapse=", ", sep="")
		return(list(valid=FALSE, probstring=probstring))	
	}
	
}

getrunjagsmethod <- function(method){

	methodmatch <- pmatch(tolower(method), c('rjags', 'simple', 'interruptible', 'parallel', 'snow', 'rjparallel', 'background', 'bgparallel', 'xgrid'))
	if(is.na(methodmatch)){
		stop(paste("Unsupported or ambiguous method '", method, "'; choose one of 'rjags', 'simple', 'interruptible', 'parallel', 'snow', 'rjparallel', 'background' or 'bgparallel'", sep=""), call.=FALSE)
	}else{
		method <- c('rjags', 'simple', 'interruptible', 'parallel', 'snow', 'rjparallel', 'background', 'bgparallel', 'xgrid')[methodmatch]
	}
	if(.Platform$OS.type=='unix' && (.Platform$GUI!="AQUA" & Sys.info()['user']=='nobody' && !(method %in% c('rjags','simple')))){
		warning("You may be trying to use a runjags method on Xgrid which won't work - choose either rjags or simple methods for using run.jags functions on Xgrid (or see the xgrid.jags functions for an alternative)")
	}
	
	return(method)
}


addmutated <- function(oldmcmc, mutate){
	
	if(!is.null(mutate)){
		
		args <- list()
		if(class(mutate)=='list'){
			if(length(mutate)==1){
				mutate <- mutate[[1]]
			}else{
				args <- mutate[-1]
				mutate <- mutate[[1]]
			}
		}
		
		if(is.character(mutate)){
		  mutate <- try(get(mutate),silent=TRUE)
			if(class(mutate)=='try-error') stop('Unable to find a function matching the character argument given to "mutate"',call.=FALSE)
		}
		
    if(!is.function(mutate))
		stop('The "mutate" argument must be a function (or a character string naming a function, or a list with first element containing the function)')
	if(length(formals(mutate))==0)
		stop('The "mutate" argument must be a function taking at least 1 argument (an mcmc matrix)')
	 
	  if(length(formals(mutate))<(1+length(args)))
      stop('The "mutate" argument must be a function with the first argument an mcmc matrix')
    
		newmcmcs <- as.mcmc.list(lapply(1:nchain(oldmcmc), function(chain){
			x <- oldmcmc[[chain]]
			args <- c(list(x), args)
			new <- do.call('mutate', args)
		
			if(is.list(new)){
				lengths <- sapply(new,length)
				if(!all(lengths == nrow(x)))
					stop('Object produced by the specified "mutate" function was invalid - a named list or matrix with an equal number of rows to that provided must be returned', call.=FALSE)
				names <- names(new)
				if(is.null(names) || any(names==""))
					stop('Object produced by the specified "mutate" function was invalid - a named list or matrix with unique variable names must be returned', call.=FALSE)
			
				new <- as.matrix(as.data.frame(new))
				dimnames(new) <- list(NULL, names)
			}
			if(is.null(dim(new))){
				dim(new) <- c(length(new), 1)
				dimnames(new) <- list(NULL, 'mutated_variable')
			}
			ds <- dim(new)
			if(length(ds)!=2 || ds[1]!=nrow(x))
				stop('Object produced by the specified "mutate" function was invalid - a named list or matrix with an equal number of rows to that provided must be returned', call.=FALSE)
			
			if(is.null(dimnames(new)[[2]]) || any(dimnames(new)[[2]]=="") || length(unique(dimnames(new)[[2]]))!=length(dimnames(new)[[2]]))
				stop('Object produced by the specified "mutate" function was invalid - a named list or matrix with unique variable names must be returned', call.=FALSE)

			if(any(dimnames(new)[[2]] %in% varnames(x))){
				notduplicated <- which(! dimnames(new)[[2]] %in% varnames(x))
				new <- new[,notduplicated,drop=FALSE]				
				if(chain==1) warning('Duplicated variables returned by the "mutate" function were removed',call.=FALSE)
			}
		
			# Copy iteration names:
			dimnames(new) <- list(dimnames(x)[[1]], dimnames(new)[[2]])

			new <- as.mcmc(cbind(x,new))			
			return(new)
		}))
		mcmc <- newmcmcs
	}else{
		mcmc <- oldmcmc
	}
	return(mcmc)
}


getsummaryargs <- function(oldsummaryargs, parentf, ignore=character(0), ...){

  # Need to first get rid of anything in summary.pars (oldsummaryargs) not in add.summary, then check ..., then overwrite with anything in ...
  
	newsummaryargs <- list(...)
	# BANOVA is currently using these so make them a note, then a warning in the next release
	newsummaryargs <- newsummaryargs[! names(newsummaryargs)%in%ignore]
	if('check.stochastic'%in%names(newsummaryargs))
		swcat('** Note:  The argument "check.stochastic" is deprecated and will be ignored **\n')
		#warning('The argument "check.stochastic" is deprecated and will be ignored',call.=FALSE)	
	newsummaryargs$check.stochastic <- NULL
	if('check.conv'%in%names(newsummaryargs))
		swcat('** Note:  The argument "check.conv" is deprecated and will be ignored **\n')
		#warning('The argument "check.stochastic" is deprecated and will be ignored',call.=FALSE)	
	newsummaryargs$check.conv <- NULL

	if('monitor.deviance'%in%names(newsummaryargs))
		warning('The argument "monitor.deviance" is deprecated and will be ignored - add "deviance" to the monitors argument instead',call.=FALSE)	
	newsummaryargs$monitor.deviance <- NULL
	if('monitor.pd'%in%names(newsummaryargs))
		warning('The argument "monitor.pd" is deprecated and will be ignored - add "full.pd" to the monitors argument instead',call.=FALSE)	
	newsummaryargs$monitor.pd <- NULL
	if('monitor.pd.i'%in%names(newsummaryargs))
		warning('The argument "monitor.pd.i" is deprecated and will be ignored - add "pd" to the monitors argument instead',call.=FALSE)	
	newsummaryargs$monitor.pd.i <- NULL
	if('monitor.popt'%in%names(newsummaryargs))
		warning('The argument "monitor.popt" is deprecated and will be ignored - add "popt" to the monitors argument instead',call.=FALSE)	
	newsummaryargs$monitor.popt <- NULL

	# Check the summary arguments:
	newsummaryargs <- getargs('add.summary', newsummaryargs, otherfnames=parentf, returnall=FALSE)
  
	# Start with add.summary defaults:
	summaryargs <- formals(add.summary)
  
	# Overwrite with any values in the old summary arguments
	summaryargs <- summaryargs[!names(summaryargs)%in%c('runjags.object',names(oldsummaryargs))]
	summaryargs <- c(summaryargs, oldsummaryargs)
	
	# Overwrite with any values passed through:
	summaryargs <- summaryargs[!names(summaryargs)%in%c(names(newsummaryargs))]
	summaryargs <- c(summaryargs, newsummaryargs)
		
	return(summaryargs)
}


weightedaverage <- function(dev1, dev2, iter1, iter2){
	
	stopifnot(all(dim(dev1)==dim(dev2)))
	
	# For combining DICs - if either are completely NA, return just NA rather than a matrix of:
	if(all(is.na(dev1)) || all(is.na(dev2)))
		return(NA)
	
	totaliter <- iter1+iter2
	dev1 <- dev1 * (iter1/totaliter)
	dev2 <- dev2 * (iter2/totaliter)
	
	return(dev1 + dev2)
	
}

checkadaptrequired <- function(rjags){
	# Adapt will return TRUE if it is complete, or an error if not (it tries to update 0):
	suppressWarnings(s <- try(adaptcomplete <- rjags::adapt(rjags, n.iter=0), silent=TRUE))
	if(class(s)=='try-error')
		adaptcomplete <- FALSE
	
	return(!adaptcomplete)
}

getstoptexts <- function(adaptfail=runjags.getOption('adapt.incomplete')%in%c('error')){
	
	st <- c('Deleting model','Adaptation incomplete','syntax error')
	
	# JAGS version 3 will stop at Adaptation incomplete with batch.jags TRUE or FALSE
	if(adaptfail || testjags(silent=TRUE)$JAGS.major==3)
		return(st)
	else
	# JAGS version 4 will not stop at Adaptation incomplete ??with batch.jags TRUE or FALSE??
		return(st[!st %in% c('Adaptation incomplete')])
	
}

loadandcheckrjags <- function(stop=TRUE, silent=FALSE){
	
	fail <- FALSE

	if(!any(.packages(TRUE)=="rjags")){
		if(!silent)
			swcat("\nThe rjags package is not installed - either install the package from CRAN or from https://sourceforge.net/projects/mcmc-jags/files/rjags/\n")
		fail <- TRUE
	}
	
	if(!fail && !requireNamespace("rjags")){
		if(!silent)
			swcat("\nThe rjags package is installed, but could not be loaded - run the testjags() function for more detailed information\n", sep="")
		fail <- TRUE
	}
	if(!fail && packageVersion('rjags') < 3.9){
		if(!silent)
			swcat("\nPlease update the rjags package to version 3-9 or later\n", call.=FALSE)
		fail <- TRUE
	}

	if(fail && stop)
		stop("Loading the rjags package failed (diagnostics are given above this error message)", call.=FALSE)
	
	return(!fail)
}
