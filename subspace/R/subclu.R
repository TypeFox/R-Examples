#'The SubClu Algorithm for Subspace Clustering
#'
#'The SubClu Algorithm follows a bottom-up framework, in which one-dimensional
#'clusters are generated with DBSCAN and then each cluster is expanded one
#'dimension at a time into a dimension that is known to have a cluster that only
#'differs in one dimension from this cluster. This expansion is done using
#'DBSCAN with the same parameters that were used for the original DBSCAN that
#'produced the clusters.
#'
#'@param data A Matrix of input data.
#'@param epsilon size of environment parameter for DBSCAN
#'@param minSupport minimum number of points parameter for DBSCAN
#'  
#'@references Karin Kailing, Hans-Peter Kriegel and Peer Kröger
#'  \emph{Density-Connected Subspace Clustering for High-Dimensional Data}
#'@examples
#'data("subspace_dataset")
#'SubClu(subspace_dataset,epsilon=1,minSupport=5)
#'
#'@family subspace clustering algorithms
#'@export
SubClu <- function(data,epsilon=4,minSupport=4) {
  arr <- java_object_from_data(data)
  #Now that the data is in the correct format, we can call into our Java Code that will then call into the
  #actual implementation of the Algorithm
  res <- rJava::.jcall("ClusteringApplier",returnSig="[Li9/subspace/base/Cluster;",method="subclu",arr,
                       epsilon,
                       as.integer(minSupport),
                       evalArray=F)
  #We can then turn the Java Clustering Object that was returned into an R-Friendly S3-Object
  res <- r_clusters_from_java_clusters(res)
  return(res)
}