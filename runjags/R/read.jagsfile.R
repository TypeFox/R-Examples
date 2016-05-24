#' @title Write a complete JAGS model to a text file
#' @name write.jagsfile
#' @aliases write.jagsfile write.JAGSfile
#' @export

#' @description
#' Writes the JAGS model, data, initial values and monitored variables etc to a file.  The model can then be run using a call to \code{link{run.jags}} with the filename as the model argument.
#' 
#' @keywords methods

#' @return
#' Returns the filename that the model was saved to (invisibly)

#' @seealso
#' \code{\link{read.jagsfile}} and \code{\link{run.jags}} for the reverse operation

#' @examples

#' # Set up a model:
#' # y = m x + c, assuming normal observation errors for y:
#' 
#' # Simulate the data
#' X <- 1:100
#' Y <- rnorm(length(X), 2*X + 10, 1)
#' 
#' # Model in the JAGS format
#' model <- "model { 
#' for(i in 1 : N){ 
#' 	Y[i] ~ dnorm(true.y[i], precision);
#' 	true.y[i] <- (m * X[i]) + c
#' } 
#' m ~ dunif(-1000,1000)
#' c ~ dunif(-1000,1000) 
#' precision ~ dexp(1)
#' }"
#' 
#' # Data and initial values in a named list format, 
#' # with explicit control over the random number
#' # generator used for each chain (optional): 
#' data <- list(X=X, Y=Y, N=length(X))
#' inits1 <- list(m=1, c=1, precision=1,
#' .RNG.name="base::Super-Duper", .RNG.seed=1)
#' inits2 <- list(m=0.1, c=10, precision=1,
#' .RNG.name="base::Wichmann-Hill", .RNG.seed=2)
#' 
#' \dontrun{
#' # Compile the model but don't update it (sample=0):
#' compiled <- run.jags(model=model, monitor=c("m", "c", "precision"), 
#' data=data, n.chains=2, inits=list(inits1,inits2), sample=0)
#' 
#' # Save the complete model to a file:
#' filepath <- write.jagsfile(compiled, file='model.txt')
#' 
#' # And run the model from the file:
#' results <- run.jags(filepath)
#' 
#' }

#' @references 
#' Lunn D, Jackson C, Best N, Thomas A, Spiegelhalter D (2012). The BUGS book: A practical introduction to Bayesian analysis. CRC press.


#' @param runjags.object a valid (but not necessarily updated) runjags object to be saved to file.  No default.

#' @param file a filename to which the model will be written.  Note that any files already existing in that location will be overwritten with no warning (see \code{\link{new_unique}} for a way to generate unique filenames).  No default.

#' @param remove.tags should the runjags tags #data#, #inits#, #monitors#, #modules#, #factories#, #residual#, #fitted# and #response# be removed from the original model code before writing it to file?  If left in, these may create conflicts with the tags automatically added to the new file.

#' @param write.data should the data also be written to file?  If FALSE, the model may not run from the file without specifying a new source of data.

#' @param write.inits should the data also be written to file?  If FALSE, the model may not run from the file without specifying new initial values.



write.jagsfile <- function(runjags.object, file, remove.tags=TRUE, write.data=TRUE, write.inits=TRUE){
	
	# Don't check valid runjags object, just check model, data, end.state, monitors, modules, factories
	if(class(runjags.object)!='runjags' || !all(c('model','data','end.state','monitor','modules','factories','response','residual','fitted')%in%names(runjags.object)))
		stop('Invalid runjags object provided')
	
	if(missing(file))
		stop('argument "file" is missing')
	s <- try(cat('', file=file, append=FALSE))
	if(class(s)=='try-error')
		stop('Unable to create the file - check directory permissions and that the filename given is valid on your system')
		
	model <- runjags.object$model
	# May need to remove #data#, #inits#, #monitor#, #modules#, #factories# from model description:
	if(remove.tags){
		modelsplit <- strsplit(model, split='[\n\r]')[[1]]
		modelsplit <- gsub('#data#.*','',modelsplit)
		modelsplit <- gsub('#inits#.*','',modelsplit)
		modelsplit <- gsub('#monitor#.*','',modelsplit)
		modelsplit <- gsub('#modules#.*','',modelsplit)
		modelsplit <- gsub('#factories#.*','',modelsplit)
		modelsplit <- gsub('#response#.*','',modelsplit)
		modelsplit <- gsub('#residual#.*','',modelsplit)
		modelsplit <- gsub('#fitted#.*','',modelsplit)
		model <- paste(modelsplit, collapse='\n')
	}	
	
	if(all(runjags.object$end.state=='') || !write.inits){
		initsline <- ''
	}else{
		initsline <- paste('\n######################################################################################################\n######################################################################################################\n#### Initial values \n######################################################################################################\n######################################################################################################\n\n', paste(paste('inits{\n', runjags.object$end.state, '}\n\n', sep=''), collapse=''), sep='')
	}
	
	magicline <- paste('#monitor# ', paste(runjags.object$monitor, collapse=', '), '\n', sep='')
	
	runjags.object$modules <- checkmodfact(runjags.object$modules, 'module')
	if(!identical(runjags.object$modules, '') && !all(is.na(runjags.object$modules)))
		magicline <- paste(magicline, '#modules# ', paste(sapply(runjags.object$modules,function(x){
			return(paste(x[1], c('off','on')[as.logical(x[2])+1], sep=' '))
		}), collapse=', '), '\n', sep='')
	
	runjags.object$factories <- checkmodfact(runjags.object$factories, 'factory')
	if(!identical(runjags.object$factories, '') && !all(is.na(runjags.object$factories)))
		magicline <- paste(magicline, '#factories# ', paste(sapply(runjags.object$factories,function(x){
			return(paste(x[1], x[2], c('off','on')[as.logical(x[3])+1], sep=' '))
		}), collapse=', '), '\n', sep='')
	
	if(!identical(runjags.object$response, '') && !all(is.na(runjags.object$response)))
		magicline <- paste(magicline, '#response# ', paste(runjags.object$response, collapse=', '), '\n', sep='')
	
	if(!identical(runjags.object$residual, '') && !all(is.na(runjags.object$residual)))
		magicline <- paste(magicline, '#residual# ', paste(runjags.object$residual, collapse=', '), '\n', sep='')
	
	if(!identical(runjags.object$fitted, '') && !all(is.na(runjags.object$fitted)))
		magicline <- paste(magicline, '#fitted# ', paste(runjags.object$fitted, collapse=', '), '\n', sep='')
	
	if(!is.null(runjags.object$summary.pars$mutate)){
		mutate <- runjags.object$summary.pars$mutate
		if(is.character(mutate)){
			s <- try(mutate <- get(mutate))
			if(class(s)=='try-error'){
				warning('Attempting to retrieve the mutate function failed; writing the JAGS file with mutate as a character instead', call.=FALSE)
				mutate <- runjags.object$summary.pars$mutate
			}
		}
		dfm <- paste(capture.output(dump('mutate', file='', evaluate=FALSE, envir=as.environment(list(mutate=mutate)))), collapse='\n')
		magicline <- paste(magicline, '\nmutate{\n', dfm, '\n}\n', sep='')
	}
	
	if(identical(runjags.object$data, '') || !write.data)
		dataline <- ''
	else
		dataline <- paste('\n\n######################################################################################################\n######################################################################################################\n#### Data \n######################################################################################################\n######################################################################################################\n\ndata{\n', runjags.object$data, '\n}\n\n', sep='')
	
	cat('######################################################################################################\n######################################################################################################\n#### JAGS model file written by runjags version ', runjagsprivate$runjagsversion, ' on ', as.character(Sys.time()), ' \n######################################################################################################\n######################################################################################################\n\n', model, '\n', magicline, initsline, dataline,  '######################################################################################################\n######################################################################################################\n', file=file, append=TRUE, sep='')
	
	invisible(file)
	
}


write.JAGSfile <- write.jagsfile



read.jagsfile <- function(file){
  
#	if(!is.character(file) || any(class(zz)=='connection'))
	if(is.runjags(file))
    	stop('Invalid model file - this argument must be specified as a character string or a character string giving a path to a file', call.=FALSE)
  
	st <- Sys.time()
	
	exists = likelystring <- logical(length(file))

	for(i in 1:length(file)){
		likelystring[i] <- any(grepl('\n',file[i],fixed=TRUE) || grepl('\r',file[i],fixed=TRUE) || (grepl('{',file[i],fixed=TRUE) && grepl('}',file[i],fixed=TRUE)))
	
		if(!likelystring[i]){
			exists[i] <- FALSE
			for(end in c('','.txt','.bug', '.R')){
				exists[i] <- suppressWarnings(try(file.exists(paste(file[i], end, sep='')), silent=TRUE))
				if(class(exists[i])=="try-error") exists[i] <- FALSE
				if(exists[i]){
					file[i] <- paste(file[i], end, sep='')
					break
				}
			}
		}
	}
	
	if(all(!likelystring) & all(!exists)){
	  if(length(file)==1){
	    error <- paste("No model file found at the file path provided: '", file.path(getwd(),file), "'", sep='')
	  }else{
	    error <- paste("No model file found at the file paths provided: ", paste("'", file.path(getwd(),file), "'", collapse=', ', sep=''),sep='')
	  }
	  stop(error)
	}

	if(length(file)==1){
		if(exists[1]) string <- paste(readLines(file, warn=FALSE), collapse="\n") else string <- file
	
	}else{
		string <- ""
		for(i in 1:length(file)){
	
			if(likelystring[i]==FALSE & exists[i]==FALSE) warning(paste("Specified file '", file[i], "' does not exist", sep=""))
			if(exists[i]) string <- paste(string, paste(readLines(file[i], warn=FALSE), collapse="\n"), sep="\n") else string <- paste(string, file[i], sep="\n")
		}
	}

	string <- paste(string, "\n", sep="")

	# Rubbish old code that didn't use gsub and was slow:
	#find <- c("\n\t", "\n ", "\r\t", "\r ", "\n\n", "\r\r") # ",\n", ",\r", ";\n", ";\r", - removed these since we need ',\n' for some winbugs data but they are added again later on after spaces between () are removed
	#for(i in 1:length(find)){
	#	repeat{
	
	#		splits <- strsplit(string, split=find[i])
	#		string <- paste(splits[[1]], collapse="\n")
	#		if(length(splits[[1]])==1) break
	
	#	}
	#}

	# New gsub based code:
	# Remove all tabs
	string <- gsub("\t", "", string, fixed=TRUE)
	# Convert carriage returns and form feeds to \n
	string <- gsub("\f", "\n", string, fixed=TRUE)
	string <- gsub("\r", "\n", string, fixed=TRUE)
	# Remove excess white space at the start of lines:
	string <- gsub("\n[[:space:]]*", "\n", string)
	
	# Get rid of all commented out lines - leaving in the special ones we need to keep:
	tokeep <- c('#bugsdata#','#Rdata#','#modeldata#','#BUGSdata#','#Rdata#')
	changeto <- c('--bugsformat--', '--rformat--', '--modelformat--', '--bugsformat--', '--rformat--')
	nohashstring <- string
	for(k in 1:length(tokeep)){
		nohashstring <- gsub(tokeep[k], paste(changeto[k],'#',sep=''), nohashstring)  # adding the hash here clears the rest of the line
	}
	nohashstring <- gsub('#[^\n]*','',nohashstring)   # \r are converted above
#	nohashstring <- paste(lapply(strsplit(nohashstring, "[\n\r]")[[1]], function(x) gsub("#.*", "", x)), collapse="\n")  #inefficient
	for(k in 1:length(tokeep)){
		nohashstring <- gsub(changeto[k], tolower(tokeep[k]), nohashstring)
	}

	# No helpful conversion of = to <- any more (was it doing this beforehand?)
	
#   MOVED TRAILING COMMA FIX ELSEWHERE
	# Remove remaining commas and semi-colons from end of lines (allowing for as many spaces as you like between, and \n)
#	nohashstring <- gsub(",[[:space:]]*\n", "\n", nohashstring)
#	nohashstring <- gsub(";[[:space:]]*\n", "\n", nohashstring)
	# the , .Dim will be put back in later
	#   MOVED TRAILING COMMA FIX ELSEWHERE
	
	if(!grepl('model[[:space:]]*\\{',string))
		stop("No valid model was found", call.=FALSE)

	# Check that the number of { and } match:
	openbraces <- as.numeric(gregexpr('\\{', nohashstring)[[1]])
	closebraces <- as.numeric(gregexpr('\\}', nohashstring)[[1]])
	if(length(openbraces)!=length(closebraces) || all(openbraces<0) || all(closebraces<0))
		stop('Unmatched number of { and } in the specified model file - ensure all (uncommented) braces are paired correctly', call.=FALSE)

	# First extract model:
	model <- paste("model{\n", winbugs.extract.big("model", string)[[1]], "\n}\n", sep="")
	
	if(length(model)>1){
		warning("More than 1 model block was found in the file - the first model was used and other(s) ignored", call.=FALSE)
		model <- model[1]
	}

	mainmutate <- winbugs.extract.big("mutate", nohashstring, remove.list=FALSE)[[1]]
	automutate <- winbugs.extract.small('mutate', string)
	
	maindata <- winbugs.extract.big("data", nohashstring)
	maindata <- sortjagsvsbugs(maindata, data.type=TRUE)	
	model <- paste(maindata$model,model,sep='')
	maindata <- maindata$fixed

	autodata <- winbugs.extract.small("data", string)

	maininits <- winbugs.extract.big("inits", nohashstring)
	maininits <- sortjagsvsbugs(maininits, data.type=FALSE)$fixed
	
	
	############### to remove ->
	if(FALSE){
	####### Was previously in extract.big but now here as I want to see the = before converting the data.  Still always want to convert inits though.
	newstring <- maininits
	for(i in 1:length(newstring)){
		temp <- strsplit(newstring[i], "")[[1]]
	
		#  Because you can't have .Dim <- structure or variable = value:
		numbers <- which(temp=="=")
		if(length(numbers)>0){
			for(k in 1:length(numbers)){
				tstring <- character(length=10)
				for(j in 1:10){
					tstring[j] <- paste(temp[pmax((numbers[k]-3-j):(numbers[k]-j), 1)], collapse="")
				}	
				if(all(tstring!=".Dim")) temp[numbers[k]] <- "<-"
			}
		}
	
		newstring[i] <- paste(temp, collapse="")
	}
	maininits <- newstring
	######


	### Also need to check for variable=value, variable=value without new line:
	for(i in 1:length(maininits)){

	# This code will change any "," to "\n" from between "<-" and "<-" unless a "(" or ")" are present:
	# [^\\)^\\(] means any character except ) or (
	# *? and perl makes it non-greedy matching
	# Need a while loop as the same character can't be used as the end point of one search and the start of another, so variable=value, variabile=value, variable=value would only remove the first (although variable=value, variable=value\nvariable=value, variable=value would be fine)
	s <- 0
	while(s!=-1){
	greps <- gregexpr("<-([^\\)^\\(]*?),([^\\)^\\(]*?)<-", maininits[i], perl=TRUE)
	s <- greps[[1]]
	e <- (s-1)+ attr(greps[[1]], "match.length")
	if(s[1]!=-1){ # If it doesn't find anything, s=-1
	finalstring <- strsplit(maininits[i], "", fixed=TRUE)[[1]]
	for(j in 1:length(s)){
	
		tstring <- finalstring[s[j]:e[j]]
		newstring <- gsub(",", "\n", tstring)
		finalstring[s[j]:e[j]] <- newstring
	
	}
	maininits[i] <- paste(finalstring, collapse="")
	}
	}
	}
	### Then make sure any new leading white space is removed (gsub is vectorised):
	maininits <- gsub("\n[[:space:]]*", "\n", maininits)
	####

	#   MOVED TRAILING COMMA FIX ELSEWHERE
	# As a final step, replace the removed , .Dim from data:
	maininits <- gsub('\n.Dim', ',\n.Dim', maininits, fixed=TRUE)
	# And remove leading and trailing white space:
	maininits <- gsub('^ *\n', '', maininits)
	maininits <- gsub('\n *\n$', '\n', maininits)
	#   MOVED TRAILING COMMA FIX ELSEWHERE

	}
	############### -> to remove
	

	autoinits <- winbugs.extract.small("inits", string)

	####### This does need to be here - an example in read.jagsfile uses it:
	monitors <- ""
	temp <- winbugs.extract.big("monitor", string)[[1]]
	for(i in 1:length(temp)){
		tempy <- strsplit(temp[i], "")[[1]]
		tempy <- paste(tempy[tempy!=" "], collapse="")
		#tempy <- paste(strsplit(temp[i], "")[[1]][strsplit(temp[i], "")!=""][[1]], collapse="")
		for(str in c("\n", "\t", ",", ":", ";")){
			tempy <- paste(strsplit(tempy, str, fixed=TRUE)[[1]], collapse="*")
		}
		tempy <- strsplit(tempy, "*", fixed=TRUE)[[1]]
		monitors <- c(monitors, tempy[tempy!=""])
	}
	####### 

	monitors <- checkvalidmonitorname(c(monitors, winbugs.extract.small("monitor", string)))
	modules <- checkmodfact(winbugs.extract.small('modules', string), 'module')
	factories <- checkmodfact(winbugs.extract.small('factories', string), 'factory')
	response <- checkvalidmonitorname(winbugs.extract.small('response', string))
	residual <- checkvalidmonitorname(winbugs.extract.small('residual', string))
	fitted <- checkvalidmonitorname(winbugs.extract.small('fitted', string))
	monitors <- monitors[monitors!=""]

	#if(length(monitors)==0) warning("No monitor blocks or tags were found")

	#if(length(maindata)==0) warning("No data blocks or tags were found")

	model[model==''] <- NA
	maindata[maindata==''] <- NA
	autodata[autodata==''] <- NA
	maininits[maininits==''] <- NA
	autoinits[autoinits==''] <- NA
	monitors[monitors==''] <- NA
  mainmutate[mainmutate==''] <- NA	
  automutate[automutate==''] <- NA	

	if(is.null(model) || length(model)==0 || all(is.na(model)))
		model <- NA
	if(is.null(maindata) || length(maindata)==0 || all(is.na(maindata)))
		maindata <- NA
	if(is.null(autodata) || length(autodata)==0 || all(is.na(autodata)))
		autodata <- NA
	if(is.null(maininits) || length(maininits)==0 || all(is.na(maininits)))
		maininits <- NA
	if(is.null(autoinits) || length(autoinits)==0 || all(is.na(autoinits)))
		autoinits <- NA
	if(is.null(monitors) || length(monitors)==0 || all(is.na(monitors)))
		monitors <- NA
	if(is.null(modules) || length(modules)==0 || all(is.na(modules)))
		modules <- NA
	if(is.null(factories) || length(factories)==0 || all(is.na(factories)))
		factories <- NA
	if(is.null(response) || length(response)==0 || all(is.na(response)))
		response <- NA
	if(is.null(residual) || length(residual)==0 || all(is.na(residual)))
		residual <- NA
	if(is.null(fitted) || length(fitted)==0 || all(is.na(fitted)))
		fitted <- NA
	if(is.null(mainmutate) || length(mainmutate)==0 || all(is.na(mainmutate)))
		mainmutate <- NA
	if(is.null(automutate) || length(automutate)==0 || all(is.na(automutate)))
		automutate <- NA

	if(identical(modules, as.character(NA)))
		modules <- ''
	if(identical(factories, as.character(NA)))
		factories <- ''
		
	output <- list(model=model, data=maindata, autodata=autodata, inits=maininits, autoinits=autoinits, monitor=monitors, modules=modules, factories=factories, response=response, residual=residual, fitted=fitted, mutate=mainmutate, automutate=automutate)

	if(runjags.getOption('debug')>=10)
		swcat('Time taken to read model file: ', timestring(st, Sys.time()), '\n', sep='')

	return(output)

}

read.JAGSfile <- read.jagsfile

read.winbugs <- read.jagsfile
read.WinBUGS <- read.winbugs


