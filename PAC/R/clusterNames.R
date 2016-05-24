#' @importFrom dplyr %>% group_by summarise_each funs
NULL

#' Get the top features that differentiate the clusters based on the maximum cluster signal levels across all samples. The maximum signal levels of all features of each cluster are found, and then the maximum signal levels are normalized across clusters. Next, within each cluster, the ranks of the features are obtained based on the normalized signal levels. Top ranked features are used to annotate each cluster
#' 
#' @param aggregatedData Variable name of aggregated data.
#' @param Num_TopFeatures The number of top features in descending order of rank of importance.
#' @return The annotation of clusters based on top features
#' @export

clusterNames<-function(aggregatedData, Num_TopFeatures = 3){
  
  ClusterID<-NULL

  cluster_index<-c(min(aggregatedData[,1]):max(aggregatedData[,1]))
  sample_names<-unique(aggregatedData[,2])
  
  data_agg<-aggregatedData[, -2]
  data_agg$count<-1
  data_agg_intermediate<-data.frame(aggregatedData[, -2] %>% group_by(ClusterID) %>% summarise_each(funs(max)))
  data_agg<-data_agg_intermediate[order(data_agg_intermediate$ClusterID),]
  data_agg[,-c(1,ncol(data_agg))]<-data_agg[,-c(1,ncol(data_agg))]
  
  M<-data_agg[,-c(1, ncol(data_agg))]
  M_colNormalized<-scale(M, center=FALSE, scale=colSums(M))
  
  topclusterFeatures<-matrix(0,nrow(M_colNormalized), Num_TopFeatures)
  for (i in 1:nrow(M_colNormalized)){
    cluster<-M_colNormalized[i,]
    topclusterFeatures[i,]<-names(cluster[order(rank(-cluster))])[1:Num_TopFeatures] #sort in descending order
  }
  
  clusterNames<-apply(topclusterFeatures[, 1:ncol(topclusterFeatures)], 1, paste, collapse="-")
  return(clusterNames)
  
}