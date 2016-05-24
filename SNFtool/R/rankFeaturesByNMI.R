## Arguments:
## data: a list, where each item in the list is a matrix of values for each data type
## W: the target network for which the NMI is calculated against for each feature
##
## Details:
## NMI is calculated based on the clustering assignments using spectral clustering
## The number of clusters is set based on the estimateNumberOfClustersGivenGraph on the target matrix 
## using default parameters.
##
## Outputs:
## A list that contains the NMI score for each feature and their ranks from highest to lowest
## output[[1]] is the NMI score
## output[[1]][[1]] is the NMI score of first data type
## output[[1]][[1]][1] is the NMI score of the first feature of the first data type
## similarly for output[[2]]... except it is the rank instead of the score

rankFeaturesByNMI <- function(data, W) 
{  
  stopifnot(class(data) == "list")
  
  NUM_OF_DATA_TYES <- length(data)
  NMI_scores <- vector(mode="list", length=NUM_OF_DATA_TYES)
  NMI_ranks <- vector(mode="list", length=NUM_OF_DATA_TYES)
  num_of_clusters_fused <- estimateNumberOfClustersGivenGraph(W)[[1]]
  clustering_fused <- spectralClustering(W, num_of_clusters_fused)
  
  for (data_type_ind in 1:NUM_OF_DATA_TYES)
  {
    NUM_OF_FEATURES <- dim(data[[data_type_ind]])[2] 
    NMI_scores[[data_type_ind]] <- vector(mode="numeric", length=NUM_OF_FEATURES)    
    
    for (feature_ind in 1:NUM_OF_FEATURES)
    {
      affinity_matrix <- affinityMatrix(
        dist2(as.matrix(data[[data_type_ind]][, feature_ind]), as.matrix(data[[data_type_ind]][, feature_ind])))      
      clustering_single_feature <- spectralClustering(affinity_matrix, num_of_clusters_fused)
      NMI_scores[[data_type_ind]][feature_ind] <- calNMI(clustering_fused, clustering_single_feature)      
    }    
    NMI_ranks[[data_type_ind]] <- rank(-NMI_scores[[data_type_ind]], ties.method="first")
  }
  
  return(list(NMI_scores, NMI_ranks))  
}