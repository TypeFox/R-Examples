#'Synthetic Subspace Data Stream
#'
#'A Random data stream that generated data points that are clustered in several
#'subspaces.
#'@param numAtts the number of dimensions of the data stream.
#'@param numCluster the average number of clusters at any point in time.
#'@param numClusterRange amount by which the actual number of clusters can
#'  deviate from numCluster.
#'@param avgSubspaceSize the average number of dimensions in the subspace of a
#'  cluster.
#'@param avgSubspaceSizeRange the amount by which the number of dimensions can
#'  deviate from avgSubspaceSize.
#'@param kernelRadii the average radii of the clusters in the model.
#'@param kernelRadiiRange the amount by which the radii can deviate from
#'  kernelRadii.
#'@param numOverlappedCluster the number of overlapped clusters at the beginning
#'  of the stream.
#'@param overlappingDegree how close the initially overlapped clusters are
#'@param densityRange how strongly the amount of points in each cluster differs
#'  from each other. 0 means all clusters have the same size. 1 is the maximum
#'  value.
#'@param noiseLevel amount of noise
#'@param noiseInCluster can noise be placed in a cluster?
#'@param speed every speed points, the clusters move by 0.01
#'@param speedRange speed/Velocity point offset
#'@param eventFrequency events happen every eventFrequency points if at least one event is enabled and numClusterRange is set
#'@param eventMergeSplit can clusters merge or split?
#'@param eventDeleteCreate can clusters be deleted or created?
#'@param subspaceEventFrequency Subspace event frequency by each cluster movement destination.
#'@param decayHorizon decay horizon
#'@param decayThreshold decay horizon threshold
#'@param modelRandomSeed number used to seed the RNG for the model
#'@param instanceRandomSeed number used to seed the RNG for the instances
#'@export
DSD_RandomRBFSubspaceGeneratorEvents <- function(
                                                  numAtts=5,  
                                                  numCluster=5,
                                                  numClusterRange=0,
                                                  avgSubspaceSize=4,
                                                  avgSubspaceSizeRange=0,
                                                  kernelRadii=0.07,
                                                  kernelRadiiRange=0,
                                                  numOverlappedCluster=0,
                                                  overlappingDegree=0,
                                                  densityRange=0,
                                                  noiseLevel=0.1,
                                                  noiseInCluster=F,
                                                  speed=200,
                                                  speedRange=0,
                                                  eventFrequency=30000,
                                                  eventMergeSplit=F,
                                                  eventDeleteCreate=F,
                                                  subspaceEventFrequency=0,
                                                  decayHorizon=1000,
                                                  decayThreshold=0.1,
                                                  modelRandomSeed=sample.int(n = (2**31)-1,1),
                                                  instanceRandomSeed=sample.int(n = (2**31)-1,1)
                                                #  Removed for now. Might be added later together with
                                                #  evaluation functionality.
                                                #  evaluationFrequency=1000,
                                                #  subEvaluation=F,
                                                #  subEvaluationFrequency=200,
                                                ) {
  if(numAtts < (avgSubspaceSize + avgSubspaceSizeRange)) stop("Parameters can not be chosen 
                                                              in such a way that subspaces 
                                                              possibly have more dimensions 
                                                              than the original data space")
  jref <- rJava::.jcall("DataStreamBuilder",
                        returnSig="Lmoa/streams/clustering/RandomRBFSubspaceGeneratorEvents;",
                        method="buildRandomRBFSubspaceGeneratorEvents",
                        as.integer(modelRandomSeed),
                        as.integer(instanceRandomSeed),
                        as.integer(numCluster),
                        as.integer(numClusterRange),
                        as.integer(avgSubspaceSize),
                        as.integer(avgSubspaceSizeRange),
                        kernelRadii,
                        kernelRadiiRange,
                        as.integer(numOverlappedCluster),
                        overlappingDegree,
                        densityRange,
                        noiseLevel,
                        as.logical(noiseInCluster),
                        as.integer(speed),
                        as.integer(speedRange),
                        as.integer(eventFrequency),
                        as.logical(eventMergeSplit),
                        as.logical(eventDeleteCreate),
                        as.integer(subspaceEventFrequency),
                        as.integer(decayHorizon),
                        decayThreshold,
                        as.integer(1000),
                        as.logical(F),
                        as.integer(100),
                        as.integer(numAtts))
  
  
  res <- structure(list(description="A synthetic data stream that generates objects in multiple clusters, each
                        of which are not necessarily clusters in the entire data space.",
                        javaObj=jref,
                        k=numCluster,
                        d=numAtts
                        ),
                   class=c("DSD_RandomRBFSubspaceGeneratorEvents","DSD_SubspaceMOA","DSD"))
  
  return(res)
}