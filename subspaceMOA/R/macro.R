#'CLIQUE algorithm for use with DSC_ThreeStage
#'
#'An implementation of the CLIQUE algorithm that can be used with
#'\link{DSC_ThreeStage}. For more details on this algorithm, consult
#'\link[subspace]{CLIQUE}
#'@param tau the density threshold used to determine whether a hypercube is
#'  dense
#'@param xi the grid size used. E.g. a value of 10 means that the dataspace is
#'  divided into 10 regions along each dimension.
#'@examples 
#'dsc <- DSC_ThreeStage(macro=DSC_clique(),micro=DSC_subspaceCluStream())
#'dsd <- DSD_RandomRBFSubspaceGeneratorEvents()
#'update(dsc,dsd,1000)
#'@export
DSC_clique <- function(xi=10,tau=0.2) {
  jref <- rJava::.jcall("MacroClustererBuilder",
                       returnSig="Lmoa/clusterers/macrosubspace/MacroSubspaceClusterer;",
                       method="buildClique",
                       as.integer(xi),tau)
  res <- structure(list(description="",xi=xi,tau=tau,javaObj=jref),
                   class=c("DSC_clique","DSC_SubspaceMOA_macro","DSC_SubspaceMOA","DSC"))
  return(res)
}
#'ProClus algorithm for use with DSC_ThreeStage
#' 
#'An implementation of the ProClus algorithm 
#'that can be used with \link{DSC_ThreeStage}.
#'For more details on this
#'algorithm, consult \link[subspace]{ProClus}.
#'@param numOfClusters Number of Clusters to be found.
#'@param avgDimensions Average number of dimensions in which each cluster
#'  resides
#'dsc <- DSC_ThreeStage(macro=DSC_proclus(),micro=DSC_subspaceCluStream())
#'dsd <- DSD_RandomRBFSubspaceGeneratorEvents()
#'update(dsc,dsd,1000)
#'@export
DSC_proclus <- function(numOfClusters=5,avgDimensions=3){
  jref <- rJava::.jcall("MacroClustererBuilder",
                       returnSig="Lmoa/clusterers/macrosubspace/MacroSubspaceClusterer;",
                       method="buildProclus",
                       as.integer(numOfClusters),as.integer(avgDimensions))
  res <- structure(list(description="",numOfClusters=numOfClusters,avgDimensions,javaObj=jref),
                   class=c("DSC_proclus","DSC_SubspaceMOA_macro","DSC_SubspaceMOA","DSC"))
  return(res)
}
#'P3C algorithm for use with DSC_ThreeStage
#' 
#'An implementation of the P3C algorithm 
#'that can be used with \link{DSC_ThreeStage}.
#'For more details on this
#'algorithm, consult \link[subspace]{P3C}
#'@param chiSquareAlpha threshold value for the chi-square distribution that is
#'  used to determine whether an area is dense.
#'@param poissonThreshold threshold value to determine wheter two bins will be
#'  merged. Note that the value provided will be used as a negative power of 10.
#'  E.g. if a value of 20 is provided here, then the algorithm will use a
#'  threshold of 1.0*10^-20.
#'dsc <- DSC_ThreeStage(macro=DSC_p3c(),micro=DSC_subspaceCluStream())
#'dsd <- DSD_RandomRBFSubspaceGeneratorEvents()
#'update(dsc,dsd,1000)
#'@export
DSC_p3c <- function(poissonThreshold=10,chiSquareAlpha=0.001) {
  jref <- rJava::.jcall("MacroClustererBuilder",
                       returnSig="Lmoa/clusterers/macrosubspace/MacroSubspaceClusterer;",
                       method="buildP3c",
                       as.integer(poissonThreshold),chiSquareAlpha)
  res <- structure(list(description="",poissonThreshold=poissonThreshold,chiSquareAlpha=chiSquareAlpha,javaObj=jref),
                   class=c("DSC_p3c","DSC_SubspaceMOA_macro","DSC_SubspaceMOA","DSC"))
  return(res)
}
#'SubClu algorithm for use with DSC_ThreeStage
#' 
#'An implementation of the SubClu algorithm 
#'that can be used with \link{DSC_ThreeStage}.
#'For more details on this
#'algorithm, consult \link[subspace]{SubClu}
#'@param epsilon this parameter determines the size of the epsilon environment for the DBSCAN that is run as a part of this algorithm.
#'@param minSupport minimum number of points in the epsilon environment.
#'@param minOutputDim minimum dimensionality that a cluster must have to be output.
#'dsc <- DSC_ThreeStage(macro=DSC_subclu(),micro=DSC_subspaceCluStream())
#'dsd <- DSD_RandomRBFSubspaceGeneratorEvents()
#'update(dsc,dsd,1000)
#'@export
DSC_subclu <- function(epsilon=0.05,minSupport=50,minOutputDim=3){
  jref <- rJava::.jcall("MacroClustererBuilder",
                       returnSig="Lmoa/clusterers/macrosubspace/MacroSubspaceClusterer;",
                       method="buildSubclu",
                       epsilon,as.integer(minSupport),as.integer(minOutputDim))
  res <- structure(list(description="",epsilon=epsilon,minSupport=minSupport,minOutputDim=minOutputDim,javaObj=jref),
                   class=c("DSC_subclu","DSC_SubspaceMOA_macro","DSC_SubspaceMOA","DSC"))
  return(res)
}