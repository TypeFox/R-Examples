#'CluStream for use with DSC_ThreeStage
#'
#'A version of the \link[streamMOA]{DSC_CluStream} algorithm that is optimized for use with \link{DSC_ThreeStage}. 
#'Do not attempt to use this as a standalone stream clustering algorithm.
#'
#'@param timeWindow window for aging the points
#'@param maxNumKernels maximum number of microclusters to be produced
#'@param kernelRadiFactor multiplier for  the kernel radius
#'@param streamSpeed number of points processed per time unit
#'@export
DSC_subspaceCluStream <- function(timeWindow=1000,
                      maxNumKernels=200,
                      kernelRadiFactor=2,
                      streamSpeed=1) {
  javaObj <- rJava::.jcall("MicroClustererBuilder",returnSig="Lmoa/clusterers/AbstractClusterer;",
                       method="buildClustream",
                       as.integer(timeWindow),
                       as.integer(maxNumKernels),
                       as.integer(kernelRadiFactor),
                       as.integer(streamSpeed))
  res <- structure(list(description="",javaObj=javaObj),
                   class=c("DSC_subspace_clustream","DSC_SubspaceMOA_micro","DSC_SubspaceMOA","DSC"))
  return(res)
  
}
#'DenStream for use with DSC_ThreeStage
#'
#'A version of the \link[streamMOA]{DSC_DenStream} algorithm that is optimized for use with \link{DSC_ThreeStage}. 
#'Do not attempt to use this as a standalone stream clustering algorithm.
#'@param horizon range of the window
#'@param epsilon defines the epsilon neighborhood
#'@param minPoints minimal number of points a cluster has to contain
#'@param beta multiplier for mu to detect outlier micro-clusters
#'@param mu minimal number of points to form a micro-cluster
#'@param initPoints number of points used to initialize the algorithm
#'@param speed number of data points processed in one time unit
#'@export
DSC_subspaceDenStream <- function(horizon=1000,
                      epsilon=0.04,
                      minPoints=4,
                      beta=0.2,
                      mu=1,
                      initPoints=1000,
                      speed=100) {
  javaObj <- rJava::.jcall("MicroClustererBuilder",returnSig="Lmoa/clusterers/AbstractClusterer;",
                       method="buildDenstream",
                       as.integer(horizon),
                       epsilon,
                       as.integer(minPoints),
                       beta,
                       mu,
                       as.integer(initPoints),
                       as.integer(speed))
  res <- structure(list(description="",javaObj=javaObj),
                   class=c("DSC_subspace_denstream","DSC_SubspaceMOA_micro","DSC_SubspaceMOA","DSC"))
  return(res)
}