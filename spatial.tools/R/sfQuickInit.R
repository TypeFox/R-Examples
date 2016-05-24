#' Quickly initializes a parallel cluster and registers it with foreach.
#' 
#' @param cpus Number of cpus.  Will default to the max available cpus.
#' @param methods Logical. Load the methods package? (if FALSE, faster startup). Default=FALSE.
#' @param ... parameters to pass to sfInit()
#' @author Jonathan A. Greenberg
#' @details (Even more) quickly start a parallel cluster with half of available
#' cpus, parallel = TRUE, and type = "SOCK" and registers it with foreach.  
#' @examples 
#' sfQuickInit(cpus=2)
#' sfQuickStop()
#' @import parallel
#' @import doParallel
#' @export

# TODO: What to do if a cluster is already running

sfQuickInit <- function(cpus,methods=FALSE,...)
{
	if(missing("cpus"))
	{
		cpus <- floor(detectCores()/2)
	}
	
	if(!is.null(getOption("spatial.tools.current.cl"))) sfQuickStop()
	
	cl <- makeCluster(spec=cpus,type="PSOCK",methods=methods)
	
	options(spatial.tools.current.cl=cl)	
	setDefaultCluster(cl=cl)
	registerDoParallel(cl)
	return(cl)
}


