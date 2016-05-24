#'@title The CLIQUE Algorithm for Subspace Clustering
#'  
#'@description The CLIQUE Algorithm finds clusters by first dividing each
#'dimension into xi equal-width intervals and saving those intervals where the
#'density is greater than tau as clusters. Then each set of two dimensions is
#'examined: If there are two intersecting intervals in these two dimensions and
#'the density in the intersection of these intervals is greater than tau, the
#'intersection is again saved as a cluster. This is repeated for all sets of
#'three, four, five,\dots dimensions. After every step adjacent clusters are 
#'replaced by a joint cluster and in the end all of the clusters are output.
#'
#'
#'@param data A Matrix of input data.
#'@param xi Number of Intervals.
#'@param tau Density Threshold.
#'  
#'@references Rakesh Agrawal, Johannes Gehrke, Dimitrios Gunopulos, and
#'  Prabhakar Raghavan. \emph{Automatic Subspace Clustering of High Dimensional
#'  Data for Data Mining Applications}. In Proc. ACM SIGMOD, 1999.
#'  
#'@examples
#'data("subspace_dataset")
#'CLIQUE(subspace_dataset,xi=40,tau=0.06)
#'@family subspace clustering algorithms
#'@export

CLIQUE  <- function(data,xi=10,tau=0.2) {
  arr <- java_object_from_data(data)
  #Now that the data is in the correct format, we can call into our Java Code that will then call into the
  #actual implementation of the CLIQUE Algorithm
  res <- rJava::.jcall("ClusteringApplier",returnSig="[Li9/subspace/base/Cluster;",method="clique",arr,as.integer(xi),tau,evalArray=F)
  #We can then turn the Java Clustering Object that was returned into an R-Friendly S3-Object
  res <- r_clusters_from_java_clusters(res)
  return(res)
}

