#'The P3C Algorithm for Projected Clustering
#'
#'The main idea of the P3C algorithm is to use statistical distributions for the
#'task of finding clusters. To this end each dimension is first split into
#'1+log_2(nrow(data)) bins and the chi-square test is used to compute the
#'probability that the sizes of these bins are uniformly distributed. If this
#'probability is bigger than 1-\emph{ChiSquareAlpha}, nothing happens. Otherwise
#'the largest bins will be removed until this is the case. The bins that were
#'removed in this way are then used to find clusters. To this end, bins that are
#'adjacent are merged. Then clusters are formed by taking a bin from one
#'dimension and determining the probability of sharing as many points as it does
#'with each bin from another dimension. The bin is then intersected with the bin
#'from another dimension where this probability is the lowest, given that it is 
#'at least lower than 1.00E-\emph{PoissonThreshold} and this is repeated until
#'no such bin is found.
#'
#'@param data A Matrix of input data.
#'@param ChiSquareAlpha probability of not being uniformly distributed that the
#'  points in  a dimension are allowed to have without assuming that there is a
#'  cluster visible from this dimension
#'@param PoissonThreshold maximum probability for a bin in another dimension to
#'  deviate from the observed bin as much as it does that is allowed. The value
#'  used for this will be 1.00*10^{-PoissonThreshold}.
#'  
#'@examples 
#'data("subspace_dataset")
#'P3C(subspace_dataset,PoissonThreshold=3)
#'@references Gabriela Moise, Jörg Sander and Martin Ester \emph{P3C: A Robust
#'  Projected Clustering Algorithm} In Proc. 6th IEEE International Conference
#'  on Data Mining 2006
#'@family subspace clustering algorithms
#'@export
P3C  <- function(data,ChiSquareAlpha=0.005,PoissonThreshold=19) {
  arr <- java_object_from_data(data)
  #Now that the data is in the correct format, we can call into our Java Code that will then call into the
  #actual implementation of the P3C Algorithm
  res <- rJava::.jcall("ClusteringApplier",returnSig="[Li9/subspace/base/Cluster;",method="p3c",arr,
                       ChiSquareAlpha,
                       as.integer(PoissonThreshold),
                       evalArray=F)
  #We can then turn the Java Clustering Object that was returned into an R-Friendly S3-Object
  res <- r_clusters_from_java_clusters(res)
  return(res)
}