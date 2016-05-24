#' @title Mutate functions to be used with runjags summary methods
#' @name mutate.functions
#' @aliases mutate.functions contrasts.mcmc contrasts.MCMC prec2sd
#' @export

#' @description
#' Objects of class \code{\link{runjags-class}} have specialised options available for print, plot and summary.  These methods allow a mutate function to be specified which produces additional variables based on the monitor variables.  These functions are examples of valid syntax, and may be useful in their own right.

#' @details
#' The contrasts.mcmc and prec2sd functions are two common applications of the mutate argument to add.summary and \code{\link{run.jags}} and can be used as examples of the expected inputs and permitted return values.  They must take an MCMC object as input, and return an MCMC object or named list with the same length.
#' This can be used to add new variables to the posterior chains that are derived from the directly monitored variables in JAGS. This allows the variables to be summarised or extracted as part of the MCMC objects as if they had been calculated in JAGS, but without the computational or storage overheads associated with calculating them in JAGS directly.  The contrasts.mcmc and prec2sd functions are examples of valid objects (but both require an argument, so will have to be passed as e.g. mutate=list('contrasts.mcmc', 'variabletocontrast')).  See the mutate argument to \code{\link{add.summary}}.

#' @keywords methods

#' @return
#' An MCMC object.

#' @seealso
#' \code{\link{add.summary}} for an applciation of these functions.

#' @param x an object of class MCMC.

#' @param vars an optional character vector of variable names.  If supplied, only variable names in the object supplied with a partial match to anything in 'vars' will be used.  Note that regular expressions are not allowed, but the caret (^) token can be used to specify the match at the start of a variable name, and a quoted vars will be matched exactly.  Default NA meaning all variables available are returned.
NULL


#' @rdname mutate.functions
contrasts.mcmc <- function(x, vars){
  if(missing(vars) || length(vars)==0)
	  stop('A character string of length >0 must be provided for "vars"')

  matches <- matchvars(checkvalidmonitorname(vars), varnames(x))
  
  if(length(matches)==0)
	  stop('The specified variables given by vars were not found in the MCMC chains')
  if(length(matches)==1)
	  stop('Only 1 variable matching vars was found in the MCMC chains - 2 or more are requried for contrasts')
  if(length(matches)>10 && !runjags.getOption('force.summary'))
	  stop('More than 10 variables matching vars were found in the MCMC chains - trying to produce this many contrasts is not recommended')
  
  if(class(x) !='mcmc')
	  stop('An object of class mcmc must be supplied')
  
  vn <- varnames(x)[matches]
  pairs <- expand.grid(one=vn, two=vn)

  use <- apply(pairs,1,function(x) return(x[1]!=x[2]))
  pairs <- unique(t(apply(pairs[use,],1,sort)))
  
  newmcmc <- vapply(1:nrow(pairs), function(r){
	  ret <- x[,as.character(pairs[r,1])] - x[,as.character(pairs[r,2])]
	  return(ret)
  }, numeric(nrow(x)))
  
  compnames <- apply(pairs,1,paste,collapse=' vs ')
  dimnames(newmcmc) <- list(dimnames(newmcmc)[[1]], compnames)
  
  return(as.mcmc(newmcmc))
}
contrasts.MCMC <- contrasts.mcmc

#' @rdname mutate.functions
prec2sd <- function(x, vars){
  if(missing(vars) || length(vars)==0)
	  stop('A character string of length >0 must be provided for "vars"')
	
  matches <- matchvars(checkvalidmonitorname(vars), varnames(x))
  
  if(length(matches)==0)
	  stop('The specified variables given by vars were not found in the MCMC chains')
  
  if(class(x) !='mcmc')
	  stop('An object of class mcmc must be supplied')
  
  vn <- varnames(x)[matches]
	
  newmcmc <- x[,vn,drop=FALSE]
  newmcmc[] <- 1/(sqrt(newmcmc[]))
  dimnames(newmcmc) <- list(dimnames(x)[[1]], paste(vn, '.sd', sep=''))
  
  return(as.mcmc(newmcmc))
}

