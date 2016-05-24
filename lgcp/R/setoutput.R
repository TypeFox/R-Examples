##' setoutput function
##'
##' Sets output functionality for \link{lgcpPredict} via the main functions \link{dump2dir} and \link{MonteCarloAverage}. Note that it is possible for
##' the user to create their own \code{gridfunction} and \code{gridmeans} schemes.
##'
##' @param gridfunction what to do with the latent field, but default this set to nothing, but could save output to a directory, see ?dump2dir
##' @param gridmeans list of Monte Carlo averages to compute, see ?MonteCarloAverage
##' @return output parameters
##' @seealso \link{lgcpPredict}, \link{dump2dir}, \link{MonteCarloAverage}
##' @export

setoutput <- function(	gridfunction=NULL,
						gridmeans=NULL){

	return(list(gridfunction=gridfunction,
				gridmeans=gridmeans))	
						
}						
