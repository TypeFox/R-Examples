#'Control options for MCMC after ADMB fitting
#'
#'Determines the options (number of steps, save interval, etc.)  for running
#'MCMC based on the estimated mode (maximum likelihood estimate) and parameter
#'variance-covariance matrix
#'
#'See the AD Model Builder reference manual. The \code{mcrb} option (reduce
#'correlation of the Hessian when constructing the candidate distribution) and
#'the \code{mcseed} options (seed for random number generator) are not yet
#'implemented; \code{mcnoscale} above may not work properly
#'
#'@param mcmc Total number of MCMC steps
#'@param mcmc2 MCMC2 steps (see ADMB-RE manual)
#'@param mcsave Thinning interval for values saved in the PSV file.  Default is
#'\code{pmax(1,floor(mcmc/1000))}, i.e. aim to save 1000 steps
#'@param mcnoscale don't rescale step size for mcmc depending on acceptance rate
#'@param mcgrope (double) Use a candidate distribution that is a mixture of a
#'multivariate normal and a fatter-tailed distribution with a proportion
#'\code{mcmcgrope} of the fatter-tailed distribution; the ADMB manual suggests
#'values of \code{mcgrope} between 0.05 and 0.1
#'@param mcmult Multiplier for the MCMC candidate distribution
#'@param mcmcpars (character) vector of parameters to track in MCMC run.
#'\emph{At least one must be specified.} ADMB produces two kinds of output for
#'MCMC.  For any \code{sdreport} parameters it will produce a \code{hst} file
#'that contains a summary histogram; \code{mcmcpars} constructs appropriate
#'\code{sdreport} parameters in the auto-generated TPL file.  Step-by-step
#'output for all parameters (regulated by \code{mcsave}) is saved in the
#'\code{PSV} file.
#'@return Returns a list of options suitable for passing as the
#'\code{mcmc.opts} argument to \code{\link{do_admb}}
#'@note Some options (\code{mcmc2}, etc.) that can be used in AD Model Builder
#'and ADMB-RE may not be available
#'@author Ben Bolker
#' @export
#'@keywords misc
#'@examples
#'
#'mcmc.control(mcmc=2000)
#'
mcmc.control <- function(mcmc=1000,
		mcmc2=0,
		mcsave,
		mcnoscale=FALSE,
		mcgrope=FALSE,
		mcmult=1,
		mcmcpars=NULL) {
	if (missing(mcsave)) mcsave <- pmax(1,floor(mcmc/1000))
	if (mcmc2>0) {
		if (missing(mcmc)) {
			mcmc <- 0
		}
		if (mcmc>0) stop("may not specify both mcmc and mcmc2>0")
	}
	r <- list(mcsave=mcsave,mcnoscale=mcnoscale,mcgrope=mcgrope,mcmult=mcmult,mcmcpars=mcmcpars)
	if (mcmc>0) c(list(mcmc=mcmc),r) else c(list(mcmc2=mcmc2),r)
}

mcmc.args <- function(L) {
	L[["mcmcpars"]] <- NULL ## don't want to include this
	argstr <- mapply(function(n,val) {
				if (is.numeric(val)) paste("-",n," ",val,sep="") else
				if (isTRUE(val)) paste("-",val,sep="")
			},names(L),L)
	paste(unlist(argstr),collapse=" ")
}
