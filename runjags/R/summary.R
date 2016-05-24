#' @title Summary statistics and plot methods for runjags class objects
#' @name add.summary
#' @aliases add.summary summary.runjags plot.runjags
#' @export
#' @export print runjags
#' @export plot runjags

#' @description
#' Objects of class \code{\link{runjags-class}} have specialised options available for print, plot and summary.  These allow various options for controlling how the output is presented, including sub-selection of variables of interest (using partial matching).

#' @details
#' The print method is designed to display option prettily, wheras the summary method is designed to return the central table (summary statistics for each variable) as a numeric matrix that can be assigned to another variable and manipulated by the user.  If summary statistics have been pre-calculated these will be returned without re-calculation by both methods, wheras only the summary method will re-calculate summary statistics if they are not already available.
#'
#' The add.summary function returns an object of class runjags, with the new summary statistics (and plots if selected) stored internally for future use.  Note that many of the summary method options can be passed to \code{\link{run.jags}} when the model is run and will be remembered for future output, although they can be modified explicitly by subsequent calls to summary or add.summary.  If the summary statistics or plots requested are identical to those stored inside the runjags object, they will not be re-calculated.  Calculation of the mode of continuous variables is possible, but requires the suggested \code{\link[modeest]{modeest}} package.

#' @keywords methods

#' @return
#' The summary method returns a numeric matrix of summary statistics for each variable (invisibly for the print method), wheras the add.summary function returns an object of class \code{\link{runjags-class}} with the new sumamry statistics (and plots if selected) stored for future use.  

#' @seealso
#' \code{\link{runjags-class}} for details on other methods available for runjags class objects

#' @param runjags.object an object of class \code{\link{runjags-class}}.

#' @param object an object of class \code{\link{runjags-class}}.

#' @param x an object of class \code{\link{runjags-class}}.

#' @param vars an optional character vector of variable names.  If supplied, only variable names in the object supplied with a partial match to anything in 'vars' will be used.  Note that regular expressions are not allowed, but the caret (^) token can be used to specify the match at the start of a variable name, and a quoted vars will be matched exactly.  Default NA meaning all variables available are returned.

#' @param mutate either a function or a list with first element a function and remaining elements arguments to this function.  This can be used to add new variables to the posterior chains that are derived from the directly monitored variables in JAGS. This allows the variables to be summarised or extracted as part of the MCMC objects as if they had been calculated in JAGS, but without the computational or storage overheads associated with calculating them in JAGS directly.  The plot, summary and as.mcmc methods for runjags objects will automatically extract the mutated variables along with the directly monitored variables.  

#' @param psrf.target the desired cutoff for 'convergence' as determined Gelman and Rubin's convergence diagnostic (see \code{\link[coda]{gelman.diag}}).  This is somewhat arbitrary, but 1.05 is a commonly used figure.

#' @param normalise.mcmc an option test transformations of the monitored variable for improved normality, which is an assumption of the Gelman and Rubin statistic.  Setting this option to FALSE will likely cause problems with calculating the psrf for highly skewed variables.

#' @param modeest.opts arguments to be passed to the \code{\link[modeest]{mlv}} function to calculate the mode of continuous variables.  Ignored if the mode.continuous option in \code{\link{runjags.options}} is set to FALSE.

#' @param confidence a numeric vector of probabilities (between 0 and 1) on which to base confidence interval calculations.

#' @param autocorr.lags a numeric vector of integers on which to base the autocorrelation diagnostic.  See also the autocorr plot type.

#' @param custom a custom function which takes a numeric object as input and outputs a single summary statistic.  This statistic will be included with the others in the print and summary method outputs.

#' @param silent.jags option to suppress feedback text produced by the summary function when summary statistics must be recalculated.

#' @param plots option to pre-draw the plots given by plot.type to facilitate more convinient assessment of convergence after the model has finished running, at the expense of requiring a larger object to stored.  The default value uses the option given in \code{\link{runjags.options}}

#' @param plot.type a character vector of plots to produce, from 'trace', 'density', 'ecdf', 'histogram', 'autocorr', 'crosscorr', 'key' or 'all'.  These are all based on the equivalent plots from the \code{\link[lattice]{lattice}} package with some modifications.

#' @param col a vector of colours to use for the different chains.  This will be used for all plot types (where relevant), including the 'key' plot which functions to label the chain numbers of the various colours.  The default uses the standard lattice colour palatte for up to 7 chains, with a rainbow palette used for larger numbers of chains, and combined chains shown in dark grey.

#' @param summary.iters the number of iterations to thin the chains to before calculating summary statistics (including all plots except the trace plot).  Setting too high a value will cause a long delay while calculating these statistics.

#' @param trace.iters the number of iterations to thin the chains to before producing traceplots.  Setting too high a value will cause large file sizes and delays displaying the trace plots.

#' @param separate.chains option to display each plot separately for different chains (except crosscorr and key).  If FALSE, either the separate chains will be shown on the same plot (for trace, density, and ecdf) or as a single plot with combined chains (for histogram and autocorr).

#' @param trace.options a list of arguments to be passed to the underlying plot function that creates the trace plots.  A colour specification should be specified using the 'col' argument above to ensure that this is the same across plot types.

#' @param density.options a list of arguments to be passed to the underlying plot function that creates the density plots.  A colour specification should be specified using the 'col' argument above to ensure that this is the same across plot types.

#' @param histogram.options a list of arguments to be passed to the underlying plot function that creates the histogram plots.  A colour specification should be specified using the 'col' argument above to ensure that this is the same across plot types.

#' @param ecdfplot.options a list of arguments to be passed to the underlying plot function that creates the ecdf plots.  A colour specification should be specified using the 'col' argument above to ensure that this is the same across plot types.

#' @param acplot.options a list of arguments to be passed to the underlying plot function that creates the autocorr plots.  A colour specification should be specified using the 'col' argument above to ensure that this is the same across plot types.

#' @param layout the layout of the runjags plots to print, a numeric vector of length 2 stating the number of rows and columns of plots.  The default value is taken from \code{\link{runjags.options}}.

#' @param file an optional filename to which plots can be saved using \code{\link[grDevices]{pdf}}.  The default "" means produce plots in the active graphics device.

#' @param new.windows option to produce each plot (or matrix of plots) on a new graphics window rather than over-writing the previous plots.  For R interfaces where plots can be cycled through (e.g. the OS X GUI and RStudio), it is likely to be preferable to produce all plots to the same device.  The default value is taken from \code{\link{runjags.options}}, which depends on the system.

#' @param digits the number of digits to display for printed numerical output.

#' @param ... additional arguments to be passed to \code{\link[grDevices]{pdf}} for the plot.runjags method, or the default print method for the print.runjags method.
NULL



#' @rdname add.summary
add.summary <- function(runjags.object, vars=NA, mutate=NA, psrf.target = 1.05, normalise.mcmc = TRUE, modeest.opts=list(), confidence=c(0.95), autocorr.lags=c(10), custom=NULL, silent.jags=runjags.getOption('silent.jags'), plots=runjags.getOption('predraw.plots'), plot.type=c('trace','ecdf','histogram','autocorr','key','crosscorr'), col=NA, summary.iters=20000, trace.iters=1000, separate.chains=FALSE, trace.options=list(), density.options=list(), histogram.options=list(), ecdfplot.options=list(), acplot.options=list()){
	
	# We may be passed some unevaluated function arguments from parent functions using getargs so evaluate everything here:
	argnames <- names(formals(add.summary))
	for(i in 1:length(argnames)){
		success <- try(assign(argnames[i], eval(get(argnames[i]))), silent=TRUE)		
		if(inherits(success, 'try-error')){
			stop(paste("object '", strsplit(as.character(success),split="'",fixed=TRUE)[[1]][2], "' not found", sep=""), call.=FALSE)
		}
	}
	
	if(identical(mutate, NA))
		mutate <- runjags.object$summary.pars$mutate
	
	runjags.object <- checkvalidrunjagsobject(runjags.object)
	if(niter(runjags.object$mcmc)<100 && !runjags.getOption('force.summary'))
		stop("Cannot produce meaningful summary statistics with less than 100 samples")
		
	# Preserve original silent.jags:
	sj <- runjags.object$summary.pars$silent.jags
	if(is.null(sj))
		sj <- silent.jags
	
	plot.type <- c("trace","density","key","ecdf","histogram","crosscorr","autocorr","all")[pmatch(tolower(plot.type), c("trace","density","key","ecdf","histogram","crosscorr","autocorr","all"))]

	if(length(plot.type)==0 || any(is.na(plot.type))) stop("The plot.type options supplied must be in 'trace', 'density', 'key', 'ecdf', 'histogram', 'autocorr', 'crosscorr' or 'all'")
	plot.type <- unique(plot.type)
	
	if(any(plot.type=="all")){
		plot.type <- c("trace","density","key","ecdf","histogram","autocorr","crosscorr")
	}
	
	# Add mutated variables to the mcmc list:
	mcmc <- addmutated(runjags.object$mcmc, mutate)	
	monitor <- runjags.object$monitor   # This is only used to see if DIC is included
	
	if(!identical(vars, NA)){
		mcmc <- combine.mcmc(mcmc, collapse.chains=FALSE, vars=vars, add.mutate=FALSE)
	}
	
	if(summary.iters < 100) stop('Too small a value provided for summary.iters - this must be >=100')
	
	if(!silent.jags)
		swcat('Calculating summary statistics...\n')
	
	thinnedmcmc <- combine.mcmc(mcmc, collapse.chains=FALSE, return.samples=min(niter(mcmc),summary.iters))
	summariesout <- runjags.summaries(fullmcmclist=mcmc, thinnedmcmclist=thinnedmcmc, psrf.target = psrf.target, normalise.mcmc = normalise.mcmc, modeest.opts=modeest.opts, confidence=confidence, autocorr.lags=autocorr.lags, custom=custom, silent=silent.jags)
	
	if(plots && all(summariesout$nonstochastic)){
		plots <- FALSE
		if(runjags.getOption('summary.warning'))
			warning('Unable to pre-draw plots - all variables are non-stochastic', call.=FALSE)
	}
	if(plots && niter(mcmc)<2){
		plots <- FALSE
		if(runjags.getOption('summary.warning'))
			warning('Unable to pre-draw plots with less than 2 iterations', call.=FALSE)
	}
	
	if(plots){
		# Don't do thin before sending to plots, otherwise the autocorr will be wrong			
		plotsout <- runjagsplots(mcmclist=thinnedmcmc[,!summariesout$nonstochastic,drop=FALSE], psrfs=summariesout$psrf$psrf[!summariesout$nonstochastic,], discrete=summariesout$discrete[!summariesout$nonstochastic], silent=silent.jags, trace='trace'%in%plot.type, density='density'%in%plot.type, histogram='histogram'%in%plot.type, ecdf='ecdf'%in%plot.type, autocorr='autocorr'%in%plot.type, crosscorr='crosscorr'%in%plot.type, key='key'%in%plot.type, col=col, trace.iters=min(trace.iters,niter(mcmc)), separate.chains=separate.chains, trace.options=trace.options, density.options=density.options, histogram.options=histogram.options, ecdfplot.options=ecdfplot.options, acplot.options=acplot.options)
	}else{
		msg <- "Predrawn plots not available - see ?add.summary or ?plot.runjags"
		plotsout <- list(trace = msg, density = msg, histogram=msg, ecdfplot=msg, key=msg, acplot=msg, ccplot=msg)
	}
	
	runjags.object <- runjags.object[! names(runjags.object)%in%names(summariesout)]
	runjags.object <- runjags.object[! names(runjags.object)%in%names(plotsout)]
	runjags.object <- runjags.object[! names(runjags.object)%in%c("summary.pars","summary.available")]
	
	newobj <- c(runjags.object, summariesout, plotsout, list(summary.available=TRUE, summary.pars=list(plots=plots, plot.type=plot.type, vars=vars, mutate=mutate, psrf.target=psrf.target, normalise.mcmc=normalise.mcmc, modeest.opts=modeest.opts, confidence=confidence, autocorr.lags=autocorr.lags, custom=custom, silent.jags=sj, col=col, trace.iters=trace.iters, separate.chains=separate.chains, trace.options=trace.options, density.options=density.options, histogram.options=histogram.options, ecdfplot.options=ecdfplot.options, acplot.options=acplot.options)))
	class(newobj) <- 'runjags'
	
	return(newobj)

}

#' @rdname add.summary
summary.runjags <- function(object, ...){
	
	# Check ... doesn't contain plots
	passed <- list(...)
	if(any(names(passed)=='plots'))
		stop('Cannot specify a "plots" argument for the summary method')
	if(!identical(list(), passed) && (is.null(names(passed)) || any(names(passed)=='')))
		stop('Unnamed extra arguments are not allowed')
  
	# If no summary saved this overrides everything:
	if(identical(object$summary.pars, NULL) || !object$summary.available || !identical(passed, list())){
		redo <- TRUE
	}else{
		# If not changing anything, we can just reprint the saved summary:
		redo <- FALSE
	}
	
	# Catches the special circumstance of previous vars=NA (all) and now a sub-selection:
	if(identical(names(passed), "vars") && !is.null(object$summaries) && !is.null(object$summary.pars$vars) && object$summary.available && is.na(object$summary.pars$vars)){
		vars <- checkvalidmonitorname(passed$vars)
		selected <- matchvars(checkvalidmonitorname(vars),  varnames(as.mcmc.list(object)))
		return(object$summaries[selected,,drop=FALSE])
	}			
	
	passed$plots <- FALSE
	passed$runjags.object <- object
	if(redo)
		object <- do.call('add.summary', passed)
	
	return(object$summaries)
}


#' @rdname add.summary
plot.runjags <- function(x, plot.type=c("trace","ecdf","histogram","autocorr","crosscorr"), vars=NA, layout=runjags.getOption('plot.layout'), new.windows=runjags.getOption('new.windows'), file="", mutate=NULL, col=NA, trace.iters=NA, separate.chains=NA, trace.options=NA, density.options=NA, histogram.options=NA, ecdfplot.options=NA, acplot.options=NA, ...){
	
	passed <- list(...)
  	if(length(passed)>0){
		unused <- names(passed)
		unused <- unused[! unused %in% names(formals(pdf))]
		if(length(unused)>0)
			stop(paste('unused argument(s) ', paste(unused,collapse=', '), sep=' '))	
	}
	
  	plot.type <- c("trace","density","key","ecdf","histogram","crosscorr","autocorr","all")[pmatch(tolower(plot.type), c("trace","density","key","ecdf","histogram","crosscorr","autocorr","all"))]

	if(length(plot.type)==0 || any(is.na(plot.type))) stop("The plot.type options supplied must be in 'trace', 'density', 'key', 'ecdf', 'histogram', 'autocorr', 'crosscorr' or 'all'")
	plot.type <- unique(plot.type)
	
	if(any(plot.type=="all")){
		plot.type <- c("trace","density","key","ecdf","histogram","autocorr","crosscorr")
	}
	
	x <- checkvalidrunjagsobject(x)
	
	if(identical(separate.chains, NA)) separate.chains <- x$summary.pars$separate.chains
 	nchains <- nchain(x$mcmc)
	
	# Need to convert logical vars to the quoted varnames here to avoid problems with un-matched lengths:
	if(is.logical(vars) && !all(is.na(vars))){
		using <- vars
		allvars <- varnames(x$mcmc)
		if(length(allvars)!=length(using)){
			stop(paste("The length of the logical vector specified to 'vars' (", length(using), ") does not match the number of monitored variables (", length(allvars), ")", sep=""), call.=FALSE)
		}
		vars <- paste("'", allvars, "'", sep="")[using]		
	}
	
	
  if(x$summary.pars$plots && all(plot.type %in% x$summary.pars$plot.type) && is.null(mutate) && identical(col, NA) && identical(trace.iters, NA) && identical(separate.chains, x$summary.pars$separate.chains) && identical(trace.options, NA) && identical(density.options, NA) && identical(histogram.options, NA) && identical(ecdfplot.options, NA) && identical(acplot.options, NA) && (identical(vars, NA) || !'crosscorr'%in%plot.type)){
		# If vars != type then we have to redo crosscorr
		nvars <- sum(!x$nonstochastic)
	}else{
		
		# If we need to re-generate plots:
		if(identical(trace.options, NA)) trace.options <- x$summary.pars$trace.options
		if(identical(density.options, NA)) density.options <- x$summary.pars$density.options
		if(identical(histogram.options, NA)) histogram.options <- x$summary.pars$histogram.options
		if(identical(ecdfplot.options, NA)) ecdfplot.options <- x$summary.pars$ecdfplot.options
		if(identical(acplot.options, NA)) acplot.options <- x$summary.pars$acplot.options
		if(identical(trace.iters, NA)) trace.iters <- x$summary.pars$trace.iters
		if(identical(col, NA)) col <- x$summary.pars$col
		
		if(!x$summary.available || !is.null(mutate)){
			# If we need to generate summaries as well (for stochastic):
			swcat('Generating summary statistics and plots (these will NOT be saved for reuse)...\n')
			x <- add.summary(x, vars=vars, mutate=mutate, psrf.target = 1.05, normalise.mcmc = FALSE, modeest.opts=list(), confidence=c(0.95), autocorr.lags=c(1), custom=NULL, silent.jags=runjags.getOption('silent.jags'), plots=TRUE, col=col, summary.iters=trace.iters, trace.iters=trace.iters, separate.chains=separate.chains, trace.options=trace.options, density.options=density.options, histogram.options=histogram.options, ecdfplot.options=ecdfplot.options, acplot.options=acplot.options)
			nvars <- sum(!x$nonstochastic)	
		}else{
			# Otherwise just get plots:

			if(all(x$nonstochastic))
				stop('Unable to produce plots - all variables are non-stochastic', call.=FALSE)
			if(niter(x$mcmc)==1)
				stop('Unable to produce plots with less than 2 iterations', call.=FALSE)
			
			selected <- matchvars(checkvalidmonitorname(vars),  names(x$nonstochastic))
			toselect <- logical(length(x$nonstochastic))
			toselect[selected] <- TRUE
			
			if(sum(toselect & !x$nonstochastic)==0)
				stop('Unable to produce plots - all selected variables are non-stochastic', call.=FALSE)
			
			# We need to re-add the mutated variables to the mcmc list - this will match the summary stats but isn't saved:
			mcmc <- addmutated(x$mcmc, x$summary.pars$mutate)
			
			nvars <- sum(!x$nonstochastic & toselect)
			if(nvars==0)
				stop('An unexpected error occured while producing plots (no plottable variables available) - please file a bug report to the package author', call.=FALSE)
			
			swcat('Generating plots...\n')
			x <- runjagsplots(mcmc[,!x$nonstochastic & toselect,drop=FALSE], psrfs=x$psrf$psrf[!x$nonstochastic & toselect,], discrete=x$discrete[!x$nonstochastic & toselect], silent=FALSE, trace='trace'%in%plot.type, density='density'%in%plot.type, histogram='histogram'%in%plot.type, ecdf='ecdf'%in%plot.type, autocorr='autocorr'%in%plot.type, crosscorr='crosscorr'%in%plot.type, key='key'%in%plot.type, col=col, trace.iters=min(trace.iters,niter(mcmc)), separate.chains=separate.chains, trace.options=trace.options, density.options=density.options, histogram.options=histogram.options, ecdfplot.options=ecdfplot.options, acplot.options=acplot.options)
		}
	}
	
	# I might update crosscorr to do separate chains at some point - for now it is always 1 (char if not done):
	stopifnot(length(x$ccplot)==1)
	nplotseach <- max(nchains * separate.chains, 1)
	
	indplots <- sum(plot.type%in%c("trace","density","autocorr","ecdf","histogram"))
	totalindplots <- indplots * nplotseach * nvars
	combplots <- ("key" %in% plot.type) + (("crosscorr" %in% plot.type) && nvars>1)

	if(nvars==1 && indplots==0 && combplots==0)
		stop('Cannot produce a crosscorr plot for 1 variable!', call.=FALSE)
		
	# Have to unlist here because the $trace is always a list with length equal to number of chains if separate.chains or 1 otherwise
	toplot <- vector('list', length=totalindplots+combplots)
	varnames=plotnames <- character(totalindplots+combplots)
	keyorcrosscorr <- rep(FALSE,(totalindplots+combplots))
	chainnames <- gsub(' ', '0', format(1:nchains, scientific=FALSE))
	
	if(separate.chains){
		chainnames <- paste('.chain_', gsub(' ', '0', format(1:nchains, scientific=FALSE)), sep='')
	}else{
		chainnames <- ''
	} 
	
	start <- 1
	if('trace'%in%plot.type){
		indexes <- start:((start-1)+(nvars*nplotseach))
		varnames[indexes] <- rep(names(x$trace), each=nplotseach)
		plotnames[indexes] <- paste(rep(names(x$trace), each=nplotseach), '.plot1', rep(chainnames, by=nvars), sep='')
		toplot[indexes] <- unlist(x$trace, recursive=FALSE)
		start <- start+(nvars*nplotseach)
	}
	if('density'%in%plot.type){
		indexes <- start:((start-1)+(nvars*nplotseach))
		varnames[indexes] <- rep(names(x$density), each=nplotseach)
		plotnames[indexes] <- paste(rep(names(x$density), each=nplotseach), '.plot2', rep(chainnames, by=nvars), sep='')
		toplot[indexes] <- unlist(x$density, recursive=FALSE)
		start <- start+(nvars*nplotseach)
	}
	if('ecdf'%in%plot.type){
		indexes <- start:((start-1)+(nvars*nplotseach))
		varnames[indexes] <- rep(names(x$ecdfplot), each=nplotseach)
		plotnames[indexes] <- paste(rep(names(x$ecdfplot), each=nplotseach), '.plot3', rep(chainnames, by=nvars), sep='')
		toplot[indexes] <- unlist(x$ecdfplot, recursive=FALSE)
		start <- start+(nvars*nplotseach)
	}
	if('histogram'%in%plot.type){
		indexes <- start:((start-1)+(nvars*nplotseach))
		varnames[indexes] <- rep(names(x$histogram), each=nplotseach)
		plotnames[indexes] <- paste(rep(names(x$histogram), each=nplotseach), '.plot4', rep(chainnames, by=nvars), sep='')
		toplot[indexes] <- unlist(x$histogram, recursive=FALSE)
		start <- start+(nvars*nplotseach)
	}
	if('autocorr'%in%plot.type){
		indexes <- start:((start-1)+(nvars*nplotseach))
		varnames[indexes] <- rep(names(x$acplot), each=nplotseach)
		plotnames[indexes] <- paste(rep(names(x$acplot), each=nplotseach), '.plot5', rep(chainnames, by=nvars), sep='')
		toplot[indexes] <- unlist(x$acplot, recursive=FALSE)
		start <- start+(nvars*nplotseach)
	}
	
	if('key'%in%plot.type){
		varnames[start] <- 'key'
		plotnames[start] <- 'key'
		toplot[start] <- x$key
		keyorcrosscorr[start] <- TRUE
		start <- start+1
	}

  	if('crosscorr'%in%plot.type && nvars>1){
		varnames[start] <- 'crosscorr'
		plotnames[start] <- 'crosscorr'
		toplot[start] <- x$ccplot
		keyorcrosscorr[start] <- TRUE
		start <- start+1
	}	
	
	if(start > 1 && indplots>0){
		# Need to sort here or the indexing gets out of whack:
    	selected <- sort(matchvars(checkvalidmonitorname(vars),  varnames[!keyorcrosscorr]))
		nvars <- length(selected)/(nplotseach*indplots)
		selected <- c(selected,which(keyorcrosscorr))
    	toplot <- toplot[selected]
		varnames <- varnames[selected]
		plotnames <- plotnames[selected]
	}else{
		selected <- (1:length(toplot))[!keyorcrosscorr]
	}
	
	names(toplot) <- plotnames
	
	# Interleave plots for the same variable:
	indexes <- numeric()
	if(length(selected)>0){
		ta <- array(1:(length(toplot)-sum(keyorcrosscorr)), dim=c(nplotseach,nvars,indplots))
  		indexes <- as.numeric(aperm(ta,c(1,3,2)))
	}
	indexes <- c(indexes, which(keyorcrosscorr))   # Will always be at the end
	toplot <- toplot[indexes]

	invisible(print.runjagsplots(toplot,layout=layout,new.windows=new.windows,file=file,...))

}

#' @rdname add.summary
#' @method print runjags
print.runjags <- function(x, vars=NA, digits = 5, ...){
	
	x <- checkvalidrunjagsobject(x)
	
	numbers <- NULL
	success <- try({
		if(!is.list(x$method) && is.na(x$method)){
			cat("\nJAGS model summary:  Model not yet updated\n\n")
		}else{
			
			chainstring <- paste("(", if(x$thin>1) "thin = ", if(x$thin>1) x$thin, if(x$thin>1) "; ", if(nchain(x$mcmc)>1) "chains = ", if(nchain(x$mcmc)>1) nchain(x$mcmc), if(nchain(x$mcmc)>1) "; ", "adapt+burnin = ", x$burnin, ")",sep="")

			if(!x$summary.available){
				cat("\nJAGS model with ", niter(x$mcmc)*nchain(x$mcmc), " samples ", chainstring, "\n", sep="") # add #monitors
				cat("\nFull summary statistics have not been pre-calculated - use either the summary method or add.summary to calculate summary statistics\n\n", sep="")
			}else{
				
				# First check that we have matching variables:
				selected <- matchvars(checkvalidmonitorname(vars), dimnames(x$summaries)[[1]])
								
				cat("\nJAGS model summary statistics from ", niter(x$mcmc)*nchain(x$mcmc), " samples ", chainstring, ":\n", sep="")

				selected <- matchvars(checkvalidmonitorname(vars),  dimnames(x$summaries)[[1]])
				numbers <- x$summaries[selected,,drop=FALSE]
			
				m <- prettifytable(numbers, digits=digits, colsequal=FALSE, nastring="--", psrfcoldollar= (x$semistochastic & niter(x$mcmc)>1))

				print.noquote(m)
				cat("\n")
				
				if(any(x$semistochastic) && niter(x$mcmc)>1)
					cat("Note: parameters marked with '$' were non-stochastic in some chains - these parameters can not be assumed to have converged!\n")				
				
				if(class(x$dic)!="character"){
          if(is.na(x$dic$dic)){
					  cat("[DIC not available from the stored object]\n", sep="")
	        }else{
            	cat("Model fit assessment:\nDIC = ", format(round(x$dic$dic, digits=digits), scientific=FALSE), if(!any(is.na(x$dic$dic.chains))) paste("  (range between chains: ", format(round(min(x$dic$dic.chains), digits=digits), scientific=FALSE), " - ", format(round(max(x$dic$dic.chains), digits=digits), scientific=FALSE), ")", sep=""), "\n", sep="")
	        }
			if(is.na(x$dic$ped)){
            	cat("[PED not available from the stored object]\n", sep="")
			}else{
				cat("PED = ", format(round(x$dic$ped, digits=digits), scientific=FALSE), if(!any(is.na(x$dic$ped.chains))) paste("  (range between chains: ", format(round(min(x$dic$ped.chains), digits=digits), scientific=FALSE), " - ", format(round(max(x$dic$ped.chains), digits=digits), scientific=FALSE), ")", sep=""), "\n", sep="")
			}
					
          if(is.na(x$dic$meanpd))
            pdstring <- ''
          else
            pdstring <- paste("pD = ", format(round(x$dic$meanpd, digits=digits), scientific=FALSE), sep="")

          if(is.na(x$dic$meanpopt))
            poptstring <- ''
          else
            poptstring <- paste("pOpt = ", format(round(x$dic$meanpopt, digits=digits), scientific=FALSE), sep="")

          cat("Estimated effective number of parameters:  ", pdstring, if(pdstring!='' && poptstring!='') ', ', poptstring, "\n\n", sep="")
				}
				
			  cat("Total time taken: ", timestring(as.numeric(x$timetaken, units="secs")), "\n\n", sep="")
			}
		}
		})
	if(inherits(success, 'try-error')) stop("An unexpected error occured in the print method for runjags class")
	invisible(numbers)
}


#' @rdname add.summary
#' @method print runjagsplots
print.runjagsplots <- function(x,layout=runjags.getOption('plot.layout'),new.windows=runjags.getOption('new.windows'),file="",...){
	
	if(!length(layout)==2) stop("The layout option must be a numeric vector of length 2")
	if(!all(layout>0)) stop("All dimensions in the layout vector must be >=1")
	
	newwindows <- new.windows
	if(Sys.getenv('RSTUDIO')!="" && newwindows){
		warning('Only a single graphics device is permitted in Rstudio - setting new.windows to FALSE')
		newwindows <- FALSE
	}
	
	if(class(x)=="runjags"){
		
		x <- checkvalidrunjagsobject(x)
		
		if(!x$summary.pars$plots){
			cat(c(x$trace, x$density, x$histogram, x$ecdf, x$acplot, x$key, x$ccplot))
			invisible(c(x$trace, x$density, x$histogram, x$ecdf, x$acplot, x$key, x$ccplot))
		}
		
		# Call plot.runjags which will do the shuffling around and recall print:
		return(plot.runjags(x, vars=NA, layout=layout, new.windows=new.windows, file=file, plot.type=x$summary.pars$plot.type, col=NA, trace.iters=NA, separate.chains=NA, trace.options=NA, density.options=NA, histogram.options=NA, ecdfplot.options=NA, acplot.options=NA))
		# vars=NA and type= here combined with return if all character above should ensure they are never re-plotted
	}
	
	if(class(x)=='runjagsplots' && !all(sapply(x,class)=='trellis')){
		# If we are directly invoking the print method for $trace we may to sort out the possible list of separate chains
		# If called from plot(rjo) then we are passed a list, and the separating of separate chains has already been done
		x <- lapply(x, function(todo){
			names(todo) <- ''   # Need to remove the names or trellis gets confused
			todo <- unlist(todo, recursive=FALSE)
			class(todo) <- 'trellis'
			return(todo)
			})
	}
  
	if(class(x)=="character"){
		cat(x)
		invisible(x)
	}
	
	if(!all(sapply(x,class)=='trellis')){
		if(runjags.getOption('debug')){
      		cat('Not all list items are trellis\n')
      		browser()
    	}
    	stop('There was an unexpected problem formatting the plots ready for printing - please contact the package author', call.=FALSE)
	}
	          
	toplot <- x

#  	if(file==""){
#		if(!newwindows){
#			swcat("Producing ", length(toplot), " plots to the active graphics device\n(see ?runjagsclass for options to this S3 method)\n",sep="")
#		}else{
#			swcat("Producing ", length(toplot), " plots to separate graphics devices\n",sep="")
#		}
#	}else{
#		swcat("Producing ", length(toplot), " plots to file\n",sep="")
#	}			
	
	# Make layout consistent with ordering of dimensions of arrays:
	layout <- layout[2:1]
	
	if(file!=""){
		grDevices::pdf(file,...)
	}else{
		passed <- list(...)
		if(length(passed)!=0)
			warning('Additional arguments to plot.runjags were ignored')
	}

	N <- length(toplot)
	if(N==1)
		layout <- c(1,1)
	if(N==2)
		layout <- c(1,2)
	
	numperpage <- layout[1] * layout[2]
	numpages <- ceiling(N/numperpage)
	
	mores <- logical(N)
	mores[] <- TRUE
	mores[pmin(N,(0:numpages)*numperpage)] <- FALSE
	
	newpage <- TRUE
	p1 <- 1
	p2 <- 1
	
	for(i in 1:N){
		if(newpage && file=="" && newwindows){
			dev.new()
		}
		class(toplot[[i]]) <- "trellis"
		print(toplot[[i]], split=c(p1,p2,layout[1],layout[2]), newpage=newpage, more=mores[i])

		newpage <- FALSE
		p1 <- p1+1
		if(p1 > layout[1]){
			p1 <- 1
			p2 <- p2+1
			if(p2 > layout[2]){
				p2 <- 1
				newpage <- TRUE
			}
		}
	}
	
	if(file!="")
		dev.off()
		
	invisible(toplot)
}

#' @rdname add.summary
#' @method plot runjagsplots
plot.runjagsplots <- print.runjagsplots
