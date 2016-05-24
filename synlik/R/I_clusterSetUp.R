
# Utility function called from withing another function to set-up a "parallel" cluster

.clusterSetUp <- function(cluster, ncores, libraries = c(), toExport = c(), exportALL = FALSE, ...)
{ 
  parentEnv <- parent.frame()
  
  # Create a cluster (if necessary) and set clusterCreated to TRUE
  if( is.null(cluster) )
  { 
    cluster <- makeCluster(ncores)
    clusterCreated <- TRUE
  } else{
    ncores <- length(cluster)
    clusterCreated <- FALSE
  }
  
  # Put the vector of names of packages I want to load in the list of stuff to export
  # I assign "libraries" to the parent environment that it can be exported by .clusterExport
  if( length(libraries > 0) ){
    toExport <- c(toExport, "libraries")
    assign("libraries", libraries, parentEnv)
  }
  
  # Load stuff in the cluster
  if( length(toExport) > 0 || exportALL ) .clusterExport(cluster = cluster, envir = parentEnv, toExport = toExport, ALL = exportALL)
  
  # Load libraries on the cluster, delete the copy of "libraries" in the parent environment
  if( length(libraries) > 0 ) {
    rm("libraries", parentEnv)
    clusterEvalQ(cluster, sapply(libraries, function(libName) invisible(require(libName, quietly = TRUE, character.only=TRUE)) ) )
  }
  
  list("cluster" = cluster, "ncores" = ncores, "clusterCreated" = clusterCreated)
}