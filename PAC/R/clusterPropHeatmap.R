#' @import data.table
#' @importFrom dplyr %>% group_by summarise_each funs
#' @import NMF
NULL

#' Make heatmap illustrating the cluster proportions across multiple samples. The aggregated data are first padded to assign size of 0 to missing clusters in some samples. Next, the numbers of events in each cluster in each sample are obtained. These values are normalized across samples to find the cluster proportions by samples. The higher the cluster proportion in one sample, the more specific the cluster is to that sample.
#' 
#' @param data_agg Variable name of aggregated data.
#' @param Colv_order Variable vector specifying the order of heatmap columns
#' @return Plots cluster proportion heatmap plot and returns a matrix of the proportion values.
#' @export

clusterPropHeatmap<-function(data_agg, Colv_order){
  
  ClusterID<-NULL
  SampleID<-NULL
  
  #pad missing clusters in aggregaged data with 0s
  cluster_index<-c(min(data_agg[,1]):max(data_agg[,1]))
  sample_names<-unique(data_agg[,2])
  cluster_index_pad<-rep(cluster_index, length(sample_names))
  sample_names_pad<-rep(sample_names, each=max(cluster_index))
  Zero_Pad<-matrix(0, length(sample_names_pad),1)
  ZeroPadding<-data.frame(cluster_index_pad,sample_names_pad, Zero_Pad)
  dat_count_samples<-data_agg[, c(1,2,ncol(data_agg))] #ID columns and cluster sizes
  colnames(ZeroPadding)<-colnames(dat_count_samples)
  padded_dat_count_samples<-rbind.data.frame(ZeroPadding,dat_count_samples)
  
  data_agg_intermediate<-padded_dat_count_samples %>% group_by(ClusterID, SampleID) %>% summarise_each(funs(sum))
  data_agg_intermediate<-data.frame(data_agg_intermediate)
  data_agg<-data_agg_intermediate[order(data_agg_intermediate$SampleID),]
  
  NumCluster<-max(data_agg[,1])
  NumSamples<-length(unique(data_agg$SampleID))  
  sample_sizes<-data_agg[,-c(1,2)]
  
  sample_total<-tapply(sample_sizes, (seq_along(sample_sizes)-1) %/% NumCluster, sum)
  sample_sizes_perc<-sample_sizes/rep(sample_total, each=NumCluster)
  
  #prepare heatmap
  diff_mat<-cbind(data_agg[,c(1,2)], sample_sizes_perc)
  colnames(diff_mat)[3]<-"clusterSizeProportion"
  dat_melt<-melt(diff_mat[1:(NumSamples*NumCluster),], id=c("SampleID", "ClusterID"))
  dat_melt_ordered<-dat_melt[order(dat_melt$SampleID), ]
  M<-matrix(dat_melt_ordered$value, nrow(dat_melt)/NumCluster, NumCluster, byrow=TRUE)
  sampleNames<-unique(dat_melt$SampleID)
  rownames(M)<-sampleNames
  
  M_colNormalized<-scale(M, center=FALSE, scale=colSums(M))
  aheatmap(M_colNormalized, Colv = Colv_order)
  
  return(M_colNormalized)
  
}