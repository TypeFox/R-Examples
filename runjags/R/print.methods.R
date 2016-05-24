#' @title Print methods for runjags helper classes
#' @name runjags.printmethods
#' @aliases runjags.printmethods print.runjagsmodel print.runjagsdata print.runjagsinits print.runjagsoutput print.mcsestats print.crosscorrstats print.gelman.with.target print.dicstats print.runjagsbginfo

#' @export

#' @description
#' Print methods for a number of classes that are associated with runjags objects, such as model, data and initial values files etc.

#' @keywords methods

#' @seealso
#' \code{\link{runjags-class}} for print and plot methods associated with the main runjags class

#' @param x the object to be printed or converted.

#' @param object the object to be summarised.

#' @param linenumbers option to display line numbers alongside model, data and initial values output (this may be helpful for debugging).  Defualt uses the option set in \code{\link{runjags.options}}.

#' @param vars an optional character vector of variable names.  If supplied, only variable names in the object supplied with a partial match to anything in 'vars' will be used.  Note that regular expressions are not allowed, but the caret (^) token can be used to specify the match at the start of a variable name, and a quoted vars will be matched exactly.  Default NA meaning all variables available are returned.

#' @param digits the number of digits to display for printed numerical output.

#' @param ... other arguments which are passed to the default print method for some methods but ignored (with/without a warning) for others
NULL



#' @rdname runjags.printmethods
#' @method print failedjags
print.failedjags <- function(x, linenumbers=runjags.getOption('linenumbers'), ...){

	if(!identical(list(), list(...)))
		warning('Additional arguments to print.failedjags are ignored')
	type <- names(x)
	if(all(sapply(x,identical,NA))){
		cat('No failed JAGS model available!\n')
	}else{
	
		if('model' %in% type && !identical(x$model, NA)){
			print.runjagsmodel(x$model, linenumbers=linenumbers)
		}
		if('data' %in% type && !identical(x$data, NA)){
			print.runjagsdata(x$data, linenumbers=linenumbers)
		}
		if('inits' %in% type && !identical(x$inits, NA)){
			print.runjagsinits(x$inits, linenumbers=linenumbers)
		}
		if('end.state' %in% type && !identical(x$end.state, NA)){
			print.runjagsinits(x$end.state, linenumbers=linenumbers)
		}
		if('output' %in% type && !identical(x$output, NA)){
			print.runjagsoutput(x$output, linenumbers=linenumbers)
		}
	}
	invisible(x)
}

#' @rdname runjags.printmethods
#' @method print runjagsmodel
print.runjagsmodel <- function(x, linenumbers=runjags.getOption('linenumbers'), ...){
	if(linenumbers){
		split <- strsplit(x,"\n",fixed=TRUE)[[1]]
		lines <- length(split)
		x <- paste(paste(format(as.character(1:lines), justify='left'),"  |  ", split, sep=""),collapse="\n")
	}
	cat(c("\nJAGS model syntax:\n\n", x, if(linenumbers) "\n\n" else "\n"),sep="")
	invisible(x)
}

#' @rdname runjags.printmethods
#' @method print runjagsdata
print.runjagsdata <- function(x, linenumbers=runjags.getOption('linenumbers'), ...){
	if(x==""){
		cat("\nNo data supplied\n")
	}else{
		if(linenumbers){
			split <- strsplit(x,"\n",fixed=TRUE)[[1]]
			lines <- length(split)
			x <- paste(paste(format(as.character(1:lines), justify='left'),"  |  ", split, sep=""),collapse="\n")
		}
		cat(c("\nJAGS data:\n\n", x, if(linenumbers) "\n\n" else "\n"),sep="")
	}
	invisible(x)
}
#' @rdname runjags.printmethods
#' @method print runjagsinits
print.runjagsinits <- function(x, linenumbers=runjags.getOption('linenumbers'), ...){
	cat(c("","JAGS chains initial values / end states:"),sep="\n")
	if(linenumbers){
		for(i in 1:length(x)){
			split <- strsplit(x[i],"\n",fixed=TRUE)[[1]]
			lines <- length(split)
			if(lines>0) x[i] <- paste(paste(format(as.character(1:lines), justify='left'),"  |  ", split, sep=""),collapse="\n")
			cat(c("\nChain ", i, ":\n\n", "", x[i],if(linenumbers) "\n"),sep="")
		}	
	}else{
		for(i in 1:length(x)){
			cat(c("\nChain ", i, ":\n", "", x[i],""),sep="")
		}
	
	}
	cat("\n")
	invisible(x)
}
#' @rdname runjags.printmethods
#' @method print runjagsoutput
print.runjagsoutput <- function(x, linenumbers=runjags.getOption('linenumbers'), ...){

	# For runjagsstudy error outputs:
	if(!is.null(names(x)) && all(grepl("simulation",names(x)))){
		cat("\nCrashed JAGS model output(s):\n")
		simnos <- as.numeric(gsub("[[:alpha:][:punct:]]","",names(x)))
		printsimno <- TRUE
		linenumbers <- FALSE		
	}else{
		cat("\nJAGS model output(s):\n")
		simnos <- 1:length(x)
		printsimno <- length(x)>1
	}
	if(linenumbers){
		for(i in 1:length(x)){
			split <- strsplit(x[i],"\n",fixed=TRUE)[[1]]
			lines <- length(split)
			if(lines>0) x[i] <- paste(paste(format(as.character(1:lines), justify='left'),"  |  ", split, sep=""),collapse="\n")
			if(printsimno) cat(c("\nSimulation ", simnos[i], ":\n\n", "", x[i],if(linenumbers) "\n"),sep="") else cat(c("\n\n", "", x[i],if(linenumbers) "\n"),sep="")
		}	
	}else{
		for(i in 1:length(x)){
			if(printsimno) cat(c("\nSimulation ", simnos[i], ":\n", "", x[i],""),sep="") else cat(c("\n", "", x[i],""),sep="")
		}
	
	}
	cat("\n")
	invisible(x)
}
#' @rdname runjags.printmethods
#' @method print rjagsoutput
print.rjagsoutput <- function(x, ...){
	# For runjagsstudy error outputs:
	if(!is.null(names(x)) && all(grepl("simulation",names(x)))){
		cat("\nCrashed rjags model output(s):\n")
		simnos <- as.numeric(gsub("[[:alpha:][:punct:]]","",names(x)))
		printsimno <- TRUE
		linenumbers <- FALSE		
	}else{
		cat("\nrjags model output(s):\n")
		simnos <- 1:length(x)
		printsimno <- length(x)>1
	}
	linenumbers <- FALSE
	if(linenumbers){
		for(i in 1:length(x)){
			split <- strsplit(x[i],"\n",fixed=TRUE)[[1]]
			lines <- length(split)
			if(lines>0) x[i] <- paste(paste(format(as.character(1:lines), justify='left'),"  |  ", split, sep=""),collapse="\n")
			if(printsimno) cat(c("\nSimulation ", simnos[i], ":\n\n", "", x[i],if(linenumbers) "\n"),sep="") else cat(c("\n\n", "", x[i],if(linenumbers) "\n"),sep="")
		}	
	}else{
		for(i in 1:length(x)){
			if(printsimno) cat(c("\nSimulation ", simnos[i], ":\n", "", x[i],""),sep="") else cat(c("\n", "", x[i],""),sep="")
		}
	
	}
	cat("\n")
	invisible(x)
}
#' @rdname runjags.printmethods
#' @method print crosscorrstats
print.crosscorrstats <- function(x, vars=NA, digits=5, ...){
    cat("Cross-correlation matrix:\n")

	selected <- matchvars(checkvalidmonitorname(vars),  dimnames(x)[[1]])
	x <- x[selected,selected,drop=FALSE]
		
	m <- prettifytable(x, digits=digits, colsequal=FALSE, nastring="")
	
	print.noquote(m)
	
	invisible(x)
}
#' @rdname runjags.printmethods
#' @method print mcsestats
print.mcsestats <- function(x, vars=NA, digits = 5, ...){
	
	x <- unlist(x)
	varnames <- names(x)
	varnames <- varnames[grepl("sseff.",varnames)]
	varnames <- gsub("sseff.","",varnames)

	selected <- matchvars(checkvalidmonitorname(vars),  varnames)
	x <- x[c(selected, selected+length(varnames), selected+(2*length(varnames)))]
	varnames <- varnames[selected]
	
    cat("Monte Carlo standard error:\n")
	numbers <- matrix(nrow=length(varnames),ncol=4,dimnames=list(varnames,c("SSeff", "SD", "MCerr", "% of SD")))	
	numbers[,1] <- as.integer(x[1:(length(x)/3)])
	numbers[,2] <- x[((length(x)/3)+1):(length(x)*2/3)]
	numbers[,3] <- x[((length(x)*2/3)+1):length(x)]
	numbers[,4] <- round(numbers[,3]/numbers[,2]*100,1)

	m <- prettifytable(numbers, digits=digits, colsequal=FALSE, nastring="")
		
	print.noquote(m)
	
    cat("\n[A rule of thumb is that the Monte Carlo error should be less than 5% of the standard deviation of the sample]\n")

	invisible(numbers)
}

#' @rdname runjags.printmethods
#' @method print gelmanwithtarget
print.gelmanwithtarget <- function(x, vars=NA, digits = 3, ...){
	
	  selected <- matchvars(checkvalidmonitorname(vars),  dimnames(x$psrf)[[1]])
	
    cat("Potential scale reduction factors:\n\n")
    print.default(x$psrf[selected,,drop=FALSE], digits = digits, ...)
    if (!is.null(x$mpsrf)) {
        cat("\nMultivariate psrf (for all monitored variables):\n\n")
        cat(format(x$mpsrf, digits = digits))
    }
    
    cat("\n\nTarget psrf\n\n")
    cat(format(x$psrf.target, digits = digits))
    cat("\n")
	invisible(x)
}
#' @rdname runjags.printmethods
#' @method print dicstats
print.dicstats <- function(x, digits=3, ...){
	if(class(x)=="character"){
		cat(x)
	}else{
		string <- paste("Deviance information criterion [mean(deviance)+mean(pd)]:  ", round(x$dic, digits=digits), "\n", sep="")
		string <- paste(string, "Penalized Expected Deviance [mean(deviance)+mean(popt)]):  ", round(x$ped, digits=digits), "\n", sep="")
			
		swcat(string)
	}
	invisible(x)
}

#' @rdname runjags.printmethods
#' @method print runjagsbginfo
print.runjagsbginfo <- function(x, ...){

	cat("\nJAGS model summary:  Model currently running in the background ... use 'results.jags(\"", x$jobname, "\")' to retrieve the results.\nStarted on ", as.character(x$startedon), " in the following directory: '", x$directory, "'\n\n", sep="")
	invisible(x)
		
}

#' @rdname runjags.printmethods
#' @method print runjagsstudy
print.runjagsstudy <- function(x,...){
	
	if(! any(c('means','singles') %in% names(x)))
		stop('Summary statistics are not available; please email the package maintainer for help!', call.=FALSE)
	
	if(identical(x$means, NA) && identical(x$singles, NA))
		stop('Summary statistics are all missing; please email the package maintainer for help!', call.=FALSE)
	
	passed <- list(...)
	if(!identical(passed, list()))
		warning('Options to the print and summary class for the runjagsstudy class are ignored', call.=FALSE)
	
	n.sims <- x$simulations-sum(x$crashed)
	dropk <- x$dropk
	if(is.null(dropk))
		dropk <- FALSE
	
	if(!identical(x$means, NA)){
		cat("\nAverage values obtained from a ", if(dropk) "drop-k" else "JAGS", " study with a ", if(n.sims==1) "single simulation" else paste("total of ", n.sims, " simulations", sep=""), if(sum(x$crashed)>0) paste(" (excluding ", sum(x$crashed), " crashed simulation", if(sum(x$crashed)>1) "s", ")", sep=""), ":\n", sep="")
		toprint <- prettifytable(x$means, colsequal=FALSE, nastring="--")
		print.noquote(toprint)
		rettab <- x$means
		if(!identical(x$singles, NA)){
			cat("\nValues obtained for variables that were stochastic for only 1 simulation:\n", sep="")
			print.noquote(prettifytable(x$singles, colsequal=FALSE),...)
			rettab <- rbind(rettab, cbind(x$singles,1))
		}		
	}else{
		if(!identical(x$singles, NA)){
			cat("\nValues obtained from a ", if(dropk) "drop-k" else "JAGS", " study with a ", if(n.sims==1) "single simulation" else paste("total of ", n.sims, " simulations", sep=""), if(sum(x$crashed)>0) paste(" (excluding ", sum(x$crashed), " crashed simulation", if(sum(x$crashed)>1) "s", ")", sep=""), ":\n", sep="")
			print.noquote(prettifytable(x$singles, colsequal=FALSE),...)
			rettab <- cbind(x$singles, Simualtions=1)
		}			
	}
	cat("\n")
	if(sum(x$crashed)>0){
		cat("The ", sum(x$crashed), " error", if(sum(x$crashed)>1) "s", " returned ", if(sum(x$crashed)>1) "have " else "has ", "been stored in the '$errors' element of the list returned from ", if(dropk) "drop.k" else "run.jags.study", "\n",sep="")
	}
	cat("Average time taken:  ", timestring(mean(as.numeric(x$timetaken, units="secs"))), " (range: ", timestring(min(as.numeric(x$timetaken, units="secs"))), " - ", timestring(max(as.numeric(x$timetaken, units="secs"))), ")\n", sep="")
	cat("Average adapt+burnin required:  ", round(mean(x$burnin)), " (range: ", round(min(x$burnin)), " - ", round(max(x$burnin)), ")\n", sep="")
	cat("Average samples required:  ", round(mean(x$sample)), " (range: ", round(min(x$sample)), " - ", round(max(x$sample)), ")\n", sep="")
	cat("\n")
	invisible(rettab)
}

#' @rdname runjags.printmethods
#' @method summary runjagsstudy
summary.runjagsstudy <- function(object,...){
	
	x <- object
	if(! any(c('means','singles') %in% names(x)))
		stop('Summary statistics are not available; please email the package maintainer for help!', call.=FALSE)
	
	if(identical(x$means, NA) && identical(x$singles, NA))
		stop('Summary statistics are all missing; please email the package maintainer for help!', call.=FALSE)
	
	passed <- list(...)
	if(!identical(passed, list()))
		warning('Options to the print and summary class for the runjagsstudy class are ignored', call.=FALSE)
	
	toret <- x$means
	if(identical(NA, toret)){
		toret <- cbind(x$singles, Simulations=1)
	}else{
		if(!identical(x$singles, NA)){
			toret <- rbind(toret, cbind(x$singles, 1))
		}
	}
	return(toret)
	
}

#' @rdname runjags.printmethods
#' @method plot runjagsstudy
plot.runjagsstudy <- function(x, ...){
	stop('There is currently no plot method assocaited with the runjagsstudy class', call.=FALSE)
}
