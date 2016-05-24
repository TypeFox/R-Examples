#' Density-based Projected Clustering over High-Dimensional Data
#' 
#' This function creates a DSC object that represents an instance of the
#' HDDStream algorithm and can be used for stream clustering.
#' @param epsilonN radius of each neighborhood
#' @param beta control the effect of mu
#' @param mu minimum number of points desired to be in a microcluster
#' @param lambda decaying parameter
#' @param initPoints number of points to use for initialization
#' @param pi number of maximal subspace dimensionality
#' @param kappa parameter to define preference weighted vector
#' @param delta defines the threshold for the variance
#' @param offline offline multiplier for epsilon
#' @param speed number of incoming points per time unit
#' 
#' @details
#' HDDStream is an algorithm for the density-based projected clustering of
#' high-dimensional data streams.
#' 
#' The algorithm is initialized by buffering the first \emph{initPoints} points that arrive and then
#' applying the \emph{PreDeCon} algorithm over these points.
#' 
#' Then, Microclusters are maintained online by adding each new point to its
#' closest core Microcluster iff doing so does not increase the projected radius of this 
#' microcluster beyond \emph{epsilonN}. If a point can not be added to a core 
#' microcluster, an attempt will be made to add it to an outlier microcluster, 
#' with the same criterion as for core microclusters. If these attempts both 
#' fail, the point will start its own microcluster. Microclusters are aged
#' according to the decaying parameter \emph{lambda}.
#' 
#' Macroclustering is performed on-demand, using the \emph{PreDeCon} algorithm.
#' @examples
#' dsc <- DSC_HDDStream()
#' dsd <- DSD_RandomRBFSubspaceGeneratorEvents()
#' update(dsc,dsd,1000)
#'@export
DSC_HDDStream <- function(epsilonN=0.1,
                          beta=0.5,
                          mu=10,
                          lambda=0.5,
                          initPoints=2000,
                          pi=30,
                          kappa=10,
                          delta=0.001,
                          offline=2,
                          speed=100) {
  javaObj <- rJava::.jcall("SubspaceClustererBuilder",
                           returnSig="Lmoa/clusterers/hddstream/HDDStream;",
                           method="buildHDDStream",
                           epsilonN,
                           beta,
                           as.integer(mu),
                           lambda,
                           as.integer(initPoints),
                           as.integer(pi),
                           as.integer(kappa),
                           delta,
                           offline,
                           as.integer(speed))
  res <- structure(list(description="HDDStream",javaObj=javaObj),
                   class=c("DSC_HDDStream","DSC_SubspaceMOA","DSC"))
  return(res)
}

