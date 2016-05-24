#'The FIRES Algorithm for Subspace Clustering
#'
#'The FIRES Algorithm follows a three phase framework: In a first phase,
#'base-clusters are generated using a clustering-algorithm on each dimension in
#'isolation. Then these base-clusters are merged in a second phase to find
#'multidimensional cluster-approximations. These approximations are then refined
#'in the third phase.
#'
#'In this implementation, the first phase consists of a run of DBSCAN with the
#'parameters \emph{base_dbscan_epsilon} and \emph{base_dbscan_minpts} on the
#'objects as they appear from each particular dimension. Then, all of these
#'base-clusters whose size is smaller than \emph{minimumpercent} of the average
#'cluster size, e.g. 25% for the standard parameter setting, are discarded
#'because they are not likely to contain important information for the
#'clustering.
#'
#'In the second phase, these base clusters are merged to produce subspace
#'cluster approximations. This is achieved by computing the
#'\emph{k}-most-similar clusters for each base-cluster. Then the set of
#'best-merge-candidates for each base-cluster is determined, which contains
#'those clusters whose \emph{k}-most-similar clusters overlap the \emph{k}-most
#'similar clusters of the cluster by at least \emph{mu}. If a cluster has at
#'least \emph{minclu} best-merge-candidates,it is considered a best-merge
#'cluster. Finally, every pair of best-merge-clusters that are
#'best-merge-candidates of each other is grouped together with all of their
#'best-merge-candidates to form the cluster approximations.
#'
#'Note that some clusters need to be split and merged with two different other
#'clusters. This is done before the merging by determining whether the
#'intersection between a cluster and its most similar cluster as well as the
#'size of the cluster without this intersection are both larger than
#'\emph{split} times the average size of the base clusters.
#'
#'In the third phase, base-clusters are repeatedly removed from
#'cluster-approximations if their removal improves the amount of objects that
#'are shared by all base-clusters in the approximation. Finally, to generate the
#'definitive clusters from the cluster approximation, for each approximation all
#'base-clusters in the approximation are combined and the a clustering algorithm
#'is performed on these points. In this implementation, DBSCAN was chosen for
#'this purpose and will be performed with the parameters 
#'\emph{post_dbscan_epsilon} and \emph{post_dbscan_minpts}.
#'
#'
#'@param data A matrix or data frame of input data.
#'@param base_dbscan_epsilon parameter for the dbscan execution that generates
#'  the base clusters
#'@param base_dbscan_minpts parameter for the dbscan execution that generates
#'  the base clusters
#'@param minimumpercent size a base-cluster must have relative to the average
#'  size of base clusters so that it is not discarded
#'@param k amount of base-clusters that every base-cluster is compared to for
#'  merging purposes
#'@param mu number of most similar clusters in which two clusters must overlap
#'  in order to be considered best-merge-clusters of each other
#'@param minclu number of best-merge-candidates a cluster must have to be
#'  considered a best-merge-cluster
#'@param split a base-cluster is split to merge it with two different clusters
#'  iff both clusters resulting from the split have at least this size in
#'  proportion to the average size of base-clusters
#'@param post_dbscan_epsilon parameter for the dbscan execution that turns the
#'  cluster-approximations into the clusters that are output at the end
#'@param post_dbscan_minpts parameter for the dbscan execution that turns the
#'  cluster-approximations into the clusters that are output at the end
#'  
#'@references Hans-Peter Kriegel, Peer Kröger, Matthias Renz and Sebastian Wurst
#'  \emph{A Generic Framework for Efficient Subspace Clustering of 
#'  High-Dimensional Data} In Proc. 5th IEEE International Conference on Data
#'  Mining, 2005
#'@examples
#'data("subspace_dataset")
#'FIRES(subspace_dataset)
#'@family subspace clustering algorithms
#'@export 
FIRES  <- function(data,base_dbscan_epsilon=1,
                        base_dbscan_minpts=4,
                        minimumpercent=25,
                        k=1,
                        mu=1,
                        minclu=1,
                        split=0.66,
                        post_dbscan_epsilon=1,
                        post_dbscan_minpts=1) {
  arr <- java_object_from_data(data)
  #Now that the data is in the correct format, we can call into our Java Code that will then call into the
  #actual implementation of the Algorithm
  res <- rJava::.jcall("ClusteringApplier",returnSig="[Li9/subspace/base/Cluster;",method="fires",arr,
                       base_dbscan_epsilon,
                       as.integer(base_dbscan_minpts),
                       minimumpercent,
                       as.integer(k),
                       as.integer(mu),
                       as.integer(minclu),
                       split,
                       post_dbscan_epsilon,
                       as.integer(post_dbscan_minpts),
                       evalArray=F)
  #We can then turn the Java Clustering Object that was returned into an R-Friendly S3-Object
  res <- r_clusters_from_java_clusters(res)
  return(res)
}
