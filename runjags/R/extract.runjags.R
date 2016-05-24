#' @title Extract peripheral information from runjags objects
#' @name extract.runjags
#' @aliases extract.runjags extract extract.jags dic.runjags

#' @export

#' @description
#' Objects of class 'runjags' are produced by \code{\link{run.jags}}, \code{\link{results.jags}} and \code{\link{autorun.jags}}, and contain the MCMC chains as well as all information required to extend the simulation.  This function allows specific information to be extracted from these functions.  For other utility methods for the runjags class, see \code{\link{runjags-class}}.

#' @details 
#' The supported options for the 'what' argument are as follows:
#' \itemize{
#' \item{crosscorr}{ - the cross-correlation matrix}
#' \item{summary}{ - the same as the summary method for runjags object}
#' \item{model}{ - the model}
#' \item{data}{ - the data}
#' \item{end.state}{ - the model state at the last iteration (or initial values for non-updated models) which will be used to start an extended simulation}
#' \item{samplers}{ - a matrix giving the sampler used for stochastic nodes (not available for all models)}
#' \item{stochastic}{ - a logical vector of length equal to the number of variables indicating which variables are stochastic, with NA values for variables that are stochastic in one chain but not others - the return value of this can be passed to the 'vars' argument for combine.mcmc etc functions}
#' \item{dic}{ - the DIC, as returned by \code{\link[rjags]{dic.samples}}}
#' \item{dic}{ - the PED, as returned by \code{\link[rjags]{dic.samples}} with type="popt"}
#' \item{sum.deviance}{ - the sum of the mean estimated deviance for each stochastic variable}
#' \item{sum.pd}{ - the sum of the mean estimated pD for each stochastic variable}
#' \item{sum.popt}{ - the sum of the mean estimated pOpt for each stochastic variable}
#' \item{mean.deviance}{ - the mean estimated pD for each stochastic variable}
#' \item{mean.pd}{ - the mean estimated pD for each stochastic variable}
#' \item{mean.popt}{ - the mean estimated pOpt for each stochastic variable}
#' \item{full.deviance}{ - the sum of the model deviance at each iteration (for each chain)}
#' \item{full.pd}{ - the sum of the estimated pD at each iteration}
#' }
#' Note that for the deviance/DIC related parameters, these will be extracted from the available information if possible, or otherwise re-sampled. 

#' @keywords models

#' @seealso
#' \code{\link{runjags-class}} for additional methods for runjags objects, \code{\link{add.summary}} for details on plot, print and summary methods for runjags class objects, \code{\link{runjags.options}} for general options available, and \code{\link{run.jags}} and \code{\link{autorun.jags}} for the functions that create objects of this class.

#' @param x an object of class runjags.

#' @param what the information contained in the runjags object to be extracted.  See the details section for the available options.

#' @param force.resample option to re-draw new deviance/DIC/PED etc samples from the model (using \code{\link[rjags]{dic.samples}}) rather than using any statistics that may already be available from the saved runjags object

#' @param ... additional options to be passed to \code{\link[rjags]{dic.samples}}
NULL
		
extract <- function(x, what, ...){
	UseMethod("extract")
}

#' @rdname extract.runjags
#' @method extract runjags
extract.runjags <- function(x, what, force.resample=FALSE, ...){
	
	if(missing(what) || length(what)!=1 || class(what)!='character')
		stop('A character vector of length 1 must be supplied for "what"')
	
	if(x$sample==0 && !what%in%c('model','data','end.state'))
		stop('The runjags object has not yet been updated - use extend.jags to update the model')
	
	what <- tolower(what)
	possibilities <- c('crosscorr', 'summary', 'model', 'data', 'end.state', 'samplers', 'stochastic', 'dic', 'ped', 'sum.deviance', 'sum.pd', 'sum.popt', 'mean.deviance', 'mean.pd', 'mean.popt', 'full.deviance', 'full.pd')
	if(!what %in% possibilities)
		stop(paste('Cannot extract "', what, '" - options are ', paste(possibilities,collapse=', '), sep=''))
	
	if(what=='crosscorr'){
	  return(x$crosscorr)
	}
	if(what=='summary'){
	  return(x$summaries)
	}
	if(what=='model'){
	  return(x$model)
	}
	if(what=='data'){
	  return(x$data)
	}
	if(what=='end.state'){
	  return(x$end.state)
	}
	if(what=='samplers'){
		if(identical(x$samplers, NA))
			stop('Samplers are not available for this model - either it does not contain any data, or all nodes were updated by forward sampling from the prior')
		else
			return(x$samplers)
	}
  	if(what=='stochastic'){
		if(!x$summary.available){
			x <- add.summary(x)
		}
  		toreturn <- x$truestochastic
		toreturn[x$semistochastic] <- NA
		return(toreturn)
  	}
  	# All the rest are DIC related:
	
	# These two are class defined by rjags, so need to do the resampling if necessary:
  	if(what=='dic'){
	  return(dic.runjags(x, type='pD', force.resample=force.resample, ...))
	}
	if(what=='ped'){
	  return(dic.runjags(x, type='popt', force.resample=force.resample, ...))
	}
	
	# Some of the rest may be short-cutable from what we have saved:
	if(what=='sum.deviance'){
		if(!force.resample && !is.na(x$deviance.sum) && !is.na(x$deviance.sum['sum.mean.deviance']))
			return(x$deviance.sum['sum.mean.deviance'])
		
		ds <- dic.runjags(x, type='pD', force.resample=force.resample, ...)
		tab <- sum(as.matrix(ds$deviance))
		names(tab) <- 'sum.mean.deviance'
		return(tab)		
	}
	
	if(what=='sum.pd'){
		if(!force.resample && !is.na(x$deviance.sum) && !is.na(x$deviance.sum['sum.mean.pD']))
			return(x$deviance.sum['sum.mean.pD'])
		
		ds <- dic.runjags(x, type='pD', force.resample=force.resample, ...)
		tab <- sum(as.matrix(ds$penalty))
		names(tab) <- 'sum.mean.pD'
		return(tab)		
	}
	
	if(what=='sum.popt'){
		if(!force.resample && !is.na(x$deviance.sum) && !is.na(x$deviance.sum['sum.mean.pOpt']))
			return(x$deviance.sum['sum.mean.pOpt'])
		
		ds <- dic.runjags(x, type='pOpt', force.resample=force.resample, ...)
		tab <- sum(as.matrix(ds$penalty))
		names(tab) <- 'sum.mean.pOpt'
		return(tab)		
	}
	
	if(what=='mean.deviance'){
		if(!force.resample && !is.na(x$deviance.table) && all(!is.na(x$deviance.table['mean.deviance'])))
			return(x$deviance.table['mean.deviance'])
		
		ds <- dic.runjags(x, type='pD', force.resample=force.resample, ...)
		tab <- t(as.matrix(ds$deviance))
		dimnames(tab) <- list('mean.deviance', dimnames(tab)[[2]])
		return(tab)		
	}
	
	if(what=='mean.pd'){
		if(!force.resample && !is.na(x$deviance.table) && all(!is.na(x$deviance.table['mean.pD'])))
			return(x$deviance.table['mean.pD'])
		
		ds <- dic.runjags(x, type='pD', force.resample=force.resample, ...)
		tab <- t(as.matrix(ds$penalty))
		dimnames(tab) <- list('mean.pD', dimnames(tab)[[2]])
		return(tab)		
	}
	
	if(what=='mean.popt'){
		if(!force.resample && !is.na(x$deviance.table) && all(!is.na(x$deviance.table['mean.pOpt'])))
			return(x$deviance.table['mean.pOpt'])
		
		ds <- dic.runjags(x, type='pOpt', force.resample=force.resample, ...)
		tab <- t(as.matrix(ds$penalty))
		dimnames(tab) <- list('mean.pOpt', dimnames(tab)[[2]])
		return(tab)		
	}
	
	if(what=='full.deviance'){
		if(any(varnames(x$mcmc)=='deviance'))
      		return(as.mcmc.list(x, vars='deviance'))
		
		stop('The full deviance can only be extracted from runjags objects originally run (or extended) with a deviance monitor')
	}
	
	if(what=='full.pd'){
		if(!identical(x$pd, NA))
			return(x$pd)
    
		stop('The full pD can only be extracted from runjags objects originally run (or extended) using the simple or interruptible methods')
	}
	
}

dic.runjags <- function(x, type='pD', force.resample=FALSE, adapt=1000, n.iter=x$sample, quiet=FALSE, ...){

	runjags.object <- checkvalidrunjagsobject(x)
	
	if(!loadandcheckrjags(FALSE))
		stop('Returning DIC objects requires the rjags package')
	
	type <- tolower(type)
	if(!type%in%c('pd','popt')) stop('The type argument must be either pD or popt')
	if(type=='pd') type <- 'pD'
	
	# If dic samples taken already, return what dic.samples would, otherwise call extend.jags
	if(force.resample || identical(runjags.object$deviance.table, NA) || n.iter!=x$sample){
		redo <- TRUE
	}else{
		redo <- FALSE
		if(tolower(type)=='pd' && any(is.na(runjags.object$deviance.table[,1:2])))
			redo <- TRUE
		if(tolower(type)=='popt' && any(is.na(runjags.object$deviance.table[,c(1,3)])))
			redo <- TRUE
		if(nrow(runjags.object$deviance.table)==1 && dimnames(runjags.object$deviance.table)[[1]]=='mean')
			redo <- TRUE
	
	}

  if(redo){
		if(!quiet) swcat('Obtaining DIC samples...\n')
		newobj <- as.jags(runjags.object, adapt=adapt, quiet=quiet)
		dics <- rjags::dic.samples(newobj, type=type, n.iter=n.iter, ...)
		if(runjags.getOption('repeatable.methods'))
			warning('The rjags object has been updated, so subsequent samples taken with extend.jags(..., method="rjags") will not match those taken using other methods!')
			
	}else{
		pencol <- if(type=='pD') 2 else 3
		dics <- list(deviance=runjags.object$deviance.table[,1], penalty=runjags.object$deviance.table[,pencol], type=type)
		class(dics) <- 'dic'
	}
	
	return(dics)
}
