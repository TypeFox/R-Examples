#' Define a parallel cluster object to be used in internal functions
#'
#' This function defines a object that is placed in a relevant internal environment defined in faoutlier.
#' Internal functions will utilize this object automatically to capitalize on parallel
#' processing architecture. The object defined is a call from \code{parallel::makeCluster()}. Note that
#' if you are defining other parallel objects (for simulation desings, for example) it is not recommended
#' to define a cluster. 
#' 
#' @aliases setCluster
#' @param ncores number of cores to be used in the returned object which is
#'   passed to \code{parallel::makeCluster()}. If no input is given the maximum number of available
#'   cores will be used
#' @param remove logical; remove previously defined cluster object?
#' @author Phil Chalmers \email{rphilip.chalmers@@gmail.com}
#' @keywords parallel
#' @export setCluster
#' @examples
#'
#' \dontrun{
#'
#' #make 4 cores available for parallel computing
#' setCluster(4)
#'
#' #' #stop and remove cores
#' setCluster(remove = TRUE)
#'
#' #use all available cores
#' setCluster()
#'
#' }
setCluster <- function(ncores, remove = FALSE){
    if(remove){
        if(is.null(faoutlierClusterEnv$CLUSTER)){
            message('There is no visible CLUSTER() definition')
            return(invisible())
        }
        parallel::stopCluster(faoutlierClusterEnv$CLUSTER)
        faoutlierClusterEnv$CLUSTER <- NULL
        faoutlierClusterEnv$ncores <- 1L
        return(invisible())
    }
    if(!is.null(faoutlierClusterEnv$CLUSTER)){
        message('CLUSTER() has already been defined')
        return(invisible())
    }
    if(missing(ncores))
        ncores <- parallel::detectCores()
    if(!is.numeric(ncores))
        stop('ncores must be numeric')
    faoutlierClusterEnv$CLUSTER <- parallel::makeCluster(ncores)
    faoutlierClusterEnv$ncores <- as.integer(ncores)
    return(invisible())
}
