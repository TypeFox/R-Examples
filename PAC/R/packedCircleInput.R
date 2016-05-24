#' @importFrom dplyr %>% group_by summarise_each funs
NULL

#' Make packed circles plot illustrating the cluster proportions across multiple samples. 
#' 
#' @param data_agg Variable name of aggregated data.
#' @param cluster_names The annotated cluster names 
#' @param CountFilter The filter threshold to remove clusters with event numbers below the threshold.
#' @param filename Set to zoomablePackedCirclesInput.txt or the name of choice for the html that draws the packed circles.
#' @return Saves the D3 packed circles plot input in the working directory.
#' @export

packedCircleInput<-function(data_agg, cluster_names, CountFilter =1000, filename="zoomablePackedCirclesInput.txt"){
  ClusterID<-NULL
  SampleID<-NULL
  
  cluster_index<-c(min(data_agg[,1]):max(data_agg[,1]))
  sample_names<-unique(data_agg[,2])
  
  cluster_index_pad<-rep(cluster_index, length(sample_names))
  sample_names_pad<-rep(sample_names, each=max(cluster_index))
  Zero_Pad<-matrix(0, length(sample_names_pad),1)
  ZeroPadding<-data.frame(cluster_index_pad,sample_names_pad, Zero_Pad)
  
  dat_count_samples<-data_agg[, c(1,2,ncol(data_agg))] #ID columns and cell clusterulation sizes
  colnames(ZeroPadding)<-colnames(dat_count_samples)
  padded_dat_count_samples<-rbind.data.frame(ZeroPadding,dat_count_samples)
  
  data_agg_intermediate<-padded_dat_count_samples %>% group_by(ClusterID, SampleID) %>% summarise_each(funs(sum))
  data_agg_intermediate<-data.frame(data_agg_intermediate)
  data_agg_padded<- data_agg_intermediate[order(data_agg_intermediate$SampleID),]
  
  cluster_names_padded<-rep(cluster_names, nrow(data_agg_padded)/length(cluster_names))
  
  toJSON_Input<-cbind(as.character(data_agg_padded[,2]), cluster_names_padded, data_agg_padded[,3])
  colnames(toJSON_Input)<-NULL  
  filteredResults_JSON(toJSON_Input, CountFilter, filename)
}
