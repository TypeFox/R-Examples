registerParallel <- function(
##title<< Set up a parallel computing front end
 pckg.parallel = 'doMC'  ##<< character string: Package to use for parallel
                         ##   computing. Has to be (for the time beeing) one of 'doMC' or 'doParallel'.
 , max.cores = 0         ##<< integer: amount of cores to use
)
  ##description<< This function automatically sets up a cluster in a consistent way
  ##              way for different parallel computing packages.
  ##details<< registerParallel sets up a cluster object of the selected package. In
  ##          principle,  this is a simple wrapper around the cluster creating
  ##          functions of these packages that provides a unified usage.
  ##seealso<< 
  ##\code{\link[foreach]{foreach}}, \code{\link[doMC]{registerDoMC}}
{
    w <- NULL
    if (max.cores == 0) 
        max.cores   <- detectCores()
    if (max.cores == 1 || getDoParWorkers() < max.cores) {
        cat(paste('Registering ', max.cores, ' cores.\n', sep = ''))
        if (pckg.parallel == 'doMC'){
            if (requireNamespace("doMC", quietly = TRUE)) {                
                w <- max.cores
                doMC::registerDoMC(w)                
            } else {
                stop(paste("Package 'doMC' does not seem to be available on this",
                           "platform. Please set argument 'pckg.parallel' to 'doParralel'"), sep = "")
            }         
        } else if (pckg.parallel == 'doParallel'){
            w <- makeCluster(max.cores)
        } else if (pckg.parallel == 'snow') {
            stop('Do not use this function to create snow clusters. Use sfInit from package snowfall instead!')
        } else {
            stop(paste('Package ', pckg.parallel, ' is not (yet) supported!', sep=''))
        }
    }
    ##value<< For 'doMC': the amount of cores and for doParallel: the cluster object created. 
    invisible(w)
}
