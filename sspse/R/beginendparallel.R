#' Parallel Processing in the \code{\link[=sspse-package]{sspse}} Package
#' 
#' As the estimation requires MCMC, \code{\link[=sspse-package]{sspse}} can take
#' advantage of multiple CPUs or CPU cores on the system on which it runs, as
#' well as computing clusters. It uses package \code{parallel} and \code{snow}
#' to facilitate this, and supports MPI cluster type and likely PSOCK.
#' 
#' The number of nodes used and the parallel API are controlled using the
#' \code{parallel} and \code{type} arguments.
#' 
#' 
#' @name beginparallel
#' @aliases beginparallel endparallel
#' @docType methods
#' @section PSOCK clusters: The \code{parallel} package is used with PSOCK
#' clusters by default, to utilize multiple cores on a system. The number of
#' cores on a system can be determined with the \code{detectCores} function.
#' 
#' This method works with the base installation of R on all platforms, and does
#' not require additional software.
#' 
#' @param parallel scale; Number of threads in which to run the sampling. 
#' Defaults to 1 (no parallelism). 
#' @param type API to use for parallel processing. Supported values are \code{"MPI"} and
#' \code{"PSOCK"}. Defaults to using the \code{parallel} package with MPI clusters.
#' @param seed integer; random number integer seed.  Defaults to \code{NULL} to
#' use whatever the state of the random number generator is at the time of the
#' @param packagenames Names of packages in which load to get the package to run
#' functions in addition to those autodetected. This argument should not be
#' needed outside of very strange setups.
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information.
#' call.
#' @examples
#' 
#' \dontrun{
#' # Uses 2 SOCK clusters for MCMLE estimation
#' N0 <- 200
#' n <- 100
#' K <- 10
#' 
#' # Create probabilities for a Waring distribution 
#' # with scaling parameter 3 and mean 5, but truncated at K=10.
#' probs <- c(0.33333333,0.19047619,0.11904762,0.07936508,0.05555556,
#'            0.04040404,0.03030303,0.02331002,0.01831502,0.01465201)
#' probs <- probs / sum(probs)
#' 
#' # Look at the degree distribution for the prior
#' # Plot these if you want
#' # plot(x=1:K,y=probs,type="l")
#' # points(x=1:K,y=probs)
#' #
#' # Create a sample
#' #
#' set.seed(1)
#' pop<-sample(1:K, size=N0, replace = TRUE, prob = probs)
#' s<-sample(pop, size=n, replace = FALSE, prob = pop)
#'  
#' out <- posteriorsize(s=s,interval=10,parallel=2)
#' plot(out, HPD.level=0.9,data=pop[s])
#' summary(out, HPD.level=0.9)
#' # Let's look at some MCMC diagnostics
#' plot(out, HPD.level=0.9,mcmc=TRUE)
#' }
#' @export
beginparallel<-function(parallel=1, type="MPI", seed=NULL, packagenames=c("sspse"),verbose=TRUE){
#   require(parallel)
    if(verbose){
     cat("Engaging warp drive using MPI ...\n")
    }
#
#   Start Cluster
#
    ### Set up the cluster
    cl <- parallel::makeCluster(parallel,type=type)
    ### initialize parallel random number streams
    if(is.null(seed)){
     parallel::clusterSetRNGStream(cl)
    }else{
     parallel::clusterSetRNGStream(cl,iseed=seed)
    }
    ### start each virtual machine with libraries loaded
    for(pkg in packagenames){
      attached <- parallel::clusterCall(cl, require, package=pkg, character.only=TRUE)      
    }
#
#   Run the jobs with Rmpi
#
    ### make sure that R has printed out console messages before go parallel
    flush.console()
    return(cl)
}
#' @export
endparallel<-function(cl, type="MPI", finalize=TRUE, verbose=TRUE){
    parallel::stopCluster(cl)
# Activate the next line for Rmpi
#   if(finalize & type=="MPI"){Rmpi::mpi.finalize()}
    invisible()
}
