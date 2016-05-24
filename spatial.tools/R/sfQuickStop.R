#' Quickly stops a parallel snowfall cluster and deregisters it from foreach.
#' 
#' @param kill Logical.  If TRUE, attempt to force-quit cluster if you get an out-of-control situation.  Default=FALSE.
#' @param ... parameters to pass to sfStop()
#' @author Jonathan A. Greenberg
#' @details (Even more) quickly stop a snowfall cluster and sets foreach back
#' to sequential mode via registerDoSEQ().
#' @examples
#' sfQuickInit(cpus=2)
#' sfQuickStop()
#' @import parallel
#' @import foreach
#' @export

sfQuickStop <- function(kill=FALSE,...)
{
	cl <- getOption("spatial.tools.current.cl")
	registerDoSEQ()
	if(kill)
	{
		pid_connections <- lapply(getOption("spatial.tools.current.cl"),function(x) return(x$con))
		sapply(pid_connections,function(x) close(x))
		options(spatial.tools.current.cl=NULL)
	} else
	{
		stopCluster(cl)
		options(spatial.tools.current.cl=NULL)
	}
	
}


