##' mcmcpars function
##'
##' A function for setting MCMC options in a run of \code{lgcpPredict} for example.
##'
##' @param mala.length default = 100,
##' @param burnin default = floor(mala.length/2),
##' @param retain thinning parameter eg operated on chain every 'retain' iteration (eg store output or compute some posterior functional)
##' @param inits optional initial values for MCMC
##' @param adaptivescheme the type of adaptive mcmc to use, see ?constanth (constant h) or ?andrieuthomsh (adaptive MCMC of Andrieu and Thoms (2008))
##' @return mcmc parameters
##' @seealso \link{lgcpPredict}
##' @export

mcmcpars <- function(	mala.length,
						burnin,
						retain,
						inits = NULL,
                        adaptivescheme){
    
    # @param MCMCdiag non-negative integer, if greater than zero saves information from the MCMC chain                    
    # if(MCMCdiag>0){
    #     warning("NOTE: MCMCdiag option is no-longer operational. For MCMC diagnostics use output.control(dump2dir='your_directory') to save chain to disk",.immediate=TRUE)
    # } 

	return(list(mala.length=mala.length,
				burnin=burnin,
				retain=retain,
				inits=inits,
                adaptivescheme=adaptivescheme))	
						
}						
