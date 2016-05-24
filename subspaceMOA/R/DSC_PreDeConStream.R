#'Density-Based Projected Clustering of Data Streams
#'
#'This function creates a DSC object that represents an instance of the
#'PreDeConStream algorithm and can be used for stream clustering.
#'
#'@param epsilonN radius of each neighborhood
#'@param beta control the effect of mu
#'@param muN minimum number of points in microclusters
#'@param muF minimum number of points in macroclusters
#'@param lambda decaying parameter
#'@param initPoints number of points to use for initialization
#'@param tau number of maximal subspace dimensionality
#'@param kappa parameter to define preference weighted vector
#'@param delta defines the threshold for the variance
#'@param offline offline multiplier for epsilon
#'@param speed processing number of incoming points per time unit
#'
#'@details 
#'The PreDeConStream algorithm is a Density-Based algorithm for the projected
#'clustering of data streams. To initially obtain a set of microclusters
#'\emph{initPoints} points are buffered and clustered using the 
#'\emph{PreDeCon} algorithm. Then, microclusters are maintained by checking for
#'each new point whether it falls within the radius of an existing microcluster,
#'similar to \link{DSC_DenStream}. Microclusters are aged according to a decay
#'paramter \emph{lambda}. Macroclusters are also maintained throughout the run
#'of the algorithm by updating the affected macroclusters, whenever a change in
#'the microcluster structure has occured, using a component of the
#'\emph{PreDeCon} algorithm to do so.
#'@examples
#'dsc <- DSC_PreDeConStream()
#'dsd <- DSD_RandomRBFSubspaceGeneratorEvents()
#'update(dsc,dsd,1000)
#'@export
DSC_PreDeConStream <- function(epsilonN=0.7,
                               beta=0.3,
                               muN=10,
                               muF=3,
                               lambda=0.1,
                               initPoints=1000,
                               tau=2,
                               kappa=10,
                               delta=0.01,
                               offline=2,
                               speed=100) {
  javaObj <- rJava::.jcall("SubspaceClustererBuilder",
                           returnSig="Lmoa/clusterers/predeconStream/PreDeConStream;",
                           method="buildPreDeConStream",
                           epsilonN,
                           beta,
                           as.integer(muN),
                           as.integer(muF),
                           lambda,
                           as.integer(initPoints),
                           as.integer(tau),
                           kappa,
                           delta,
                           offline,
                           as.integer(speed))
  res <- structure(list(description="PreDeConStream",javaObj=javaObj),
                   class=c("DSC_PreDeConStream","DSC_SubspaceMOA","DSC"))
  return(res)
}