#'@title The ProClus Algorithm for Projected Clustering
#'  
#'@description The ProClus algorithm works in a manner similar to K-Medoids. 
#'  Initially, a set of medoids of a size that is proportional to k is chosen. 
#'  Then medoids that are likely to be outliers or are part of a cluster that is
#'  better represented by another medoid are removed until k medoids are left. 
#'  Clusters are then assumed to be around these medoids.
#'  
#'@param data A Matrix of input data.
#'@param k Number of Clusters to be found.
#'@param d Average number of dimensions in which the clusters reside
#'@references C. C. Aggarwal and C. Procopiuc \emph{Fast Algorithms for 
#'  Projected Clustering}. In Proc. ACM SIGMOD 1999.
#'  
#'  
#'@examples 
#'data("subspace_dataset")
#'ProClus(subspace_dataset,k=12,d=2.5)
#'@family subspace clustering algorithms
#'@export
ProClus <- function(data,k=4,d=3) {
  arr <- java_object_from_data(data)
  #Now that the data is in the correct format, we can call into our Java Code that will then call into the
  #actual implementation of the Algorithm
  res <- rJava::.jcall("ClusteringApplier",returnSig="[Li9/subspace/base/Cluster;",method="proclus",arr,
                        as.integer(k),
                       as.integer(d),
                       evalArray=F)
  #We can then turn the Java Clustering Object that was returned into an R-Friendly S3-Object
  res <- r_clusters_from_java_clusters(res)
  return(res)
}