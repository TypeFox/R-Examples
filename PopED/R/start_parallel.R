#' Start parallel computational processes
#' 
#' This tool chooses the type of parallelization process to use based on the
#' computer OS being used.  For windows the default is "snow" and for Linux-like
#' systems the default is "multicore"
#' 
#' @param parallel Should the parallel functionality start up?
#' @param num_cores How many cores to use.  Default is
#'   \code{parallel::detectCores()}. See \code{\link[parallel]{detectCores}} for more information.
#' @param parallel_type Which type of parallelization should be used? Can be
#'   "snow" or "multicore".  "snow"  works on Linux-like systems & Windows.
#'   "multicore" works only on Linux-like systems.  By default this is chosen
#'   for you depending on your operating system.
#' @param seed The random seed to use.
#' @param dlls If the computations require compiled code (DLL's) and you are
#'   using the "snow" method then you need to specify the name of the DLL's without 
#'   the extension as a text vector \code{c("this_file","that_file")}. 
#' @param ... Arguments passed to \code{\link[parallel]{makeCluster}}
#'   
#' @inheritParams optim_LS

#'
#' @return An atomic vector (TRUE or FALSE) with two attributes: "type" and "cores".
#'
#' 
#' @export
start_parallel <- function(parallel=TRUE,
                           num_cores=NULL,
                           parallel_type=NULL,
                           seed=NULL,
                           dlls=NULL,
                           ...)
{
  # Start parallel computing for poped package
  # edited from GA package version startParallel.R
  
  if(is.null(num_cores)) num_cores <- parallel::detectCores()
  if(is.null(parallel_type)) parallel_type <- if(.Platform$OS.type == "windows") 
    "snow" else "multicore"
  attr(parallel, "type") <- parallel_type
  attr(parallel, "cores") <- num_cores
  
  # start "parallel backend" if needed
  if(parallel){ 
    if(parallel_type == "snow"){ 
      # snow functionality on Unix-like systems & Windows
      cl <- parallel::makeCluster(num_cores, ...)
      attr(parallel, "cluster") <- cl
      
      # export parent environment
      varlist <- ls(envir = parent.frame(), all.names = TRUE)
      varlist <- varlist[varlist != "..."]
      list(...) #evaluate any promises
      parallel::clusterExport(cl, varlist = varlist,
                              # envir = parent.env(environment())
                              envir = parent.frame() )
      
      # export global environment (workspace)
      parallel::clusterExport(cl, 
                              varlist = ls(envir = globalenv(), 
                                           all.names = TRUE),
                              envir = globalenv())
      
      # load current packages in workers
      pkgs <- .packages()
      foo <- lapply(pkgs, function(pkg) 
        parallel::clusterCall(cl, library, package = pkg, 
                              character.only = TRUE))
      if(!is.null(seed)) parallel::clusterSetRNGStream(cl, seed)
      
      if(!is.null(dlls)){
        for(i in dlls){
          parallel::clusterCall(cl, 
                                dyn.load,
                                x=paste0(i,.Platform$dynlib.ext))
        }
      }
      
    } else if(parallel_type == "multicore") { 
      if(!is.null(seed)){
        RNGkind("L'Ecuyer-CMRG") 
        set.seed(seed)
        #set.seed(seed, "L'Ecuyer")
      } 
      # multicore functionality on Unix-like systems
    }
    else { stop("Only 'snow' and 'multicore' clusters allowed!") }
  }
  
  return(parallel)
}