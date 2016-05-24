#'@export
get_centers.DSC_SubspaceMOA <- function(x,type=c("auto","micro","macro"),...) {
  if("auto" %in% type) type <- "macro"  
  if("macro" %in% type) {
    macro_clustering <- rJava::.jcall(x$javaObj,"Lmoa/cluster/SubspaceClustering;","getClusteringResult")
    if(rJava::is.jnull(macro_clustering)) return(NULL)
    res <- rJava::.jcall("ClusteringAccessor",
                         returnSig="[[D",method="getClusteringResult",
                         macro_clustering,simplify=T)
    if(length(as.vector(res))==1 & as.vector(res)[1]==0)return(NULL)
    return(data.frame(res))
  } else if("micro" %in% type) {
    cast <- rJava::.jcast(x$javaObj,"moa/clusterers/AbstractSubspaceClusterer")
    res <- rJava::.jcall("ClusteringAccessor",returnSig="[[D",method="getMicroClusteringResult",cast,simplify = T)
    if(length(as.vector(res))==1 & as.vector(res)[1]==0)return(NULL)
    return(data.frame(res))
  } else {
    stop("Not implemented yet")
  }
}
#'@export
get_weights.DSC_SubspaceMOA <- function(x,type=c("auto","micro","macro"),scale=NULL,...) {
  if("auto" %in% type) type <- "macro"  
  if("macro" %in% type) {
    macro_clustering <- rJava::.jcall(x$javaObj,"Lmoa/cluster/SubspaceClustering;","getClusteringResult")
    if(rJava::is.jnull(macro_clustering)) return(NULL)
    res <- rJava::.jcall("ClusteringAccessor",returnSig="[D",method="getMacroClusteringWeights",
                         macro_clustering,simplify=T)
    if(length(as.vector(res))==1 & as.vector(res)[1]==0)return(NULL)
    return(scale_weights(res,scale=scale))
  } else if("micro" %in% type) {
    cast <- rJava::.jcast(x$javaObj,"moa/clusterers/SubspaceClusterer")
    micro_clustering <- rJava::.jcall(x$javaObj,"Lmoa/cluster/Clustering;","getMicroClusteringResult")
    if(rJava::is.jnull(micro_clustering)) return(NULL)
    res <- rJava::.jcall("ClusteringAccessor",returnSig="[D",method="getMicroClusteringWeights",micro_clustering,simplify=T)
    if(length(as.vector(res))==1 & as.vector(res)[1]==0)return(NULL)
    return(scale_weights(res,scale=scale))
  } else {
    stop("Not implemented yet")
  }
}

#'@export
#'@import stream
update.DSC_SubspaceMOA <- function(object,dsd,n = 1, verbose = FALSE, ...) {
  points <- get_points(dsd,n)
  dsc <- object
  apply(points,1,function(row){
    instance <- rJava::.jnew("moa/core/SubspaceInstance",1.0,row)
    wekaInstance <- rJava::.jcast(instance,new.class="weka/core/Instance")
    rJava::.jcall(dsc$javaObj,returnSig="V",method="trainOnInstance",wekaInstance)
  })
}
