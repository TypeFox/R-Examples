#' @importFrom dplyr %>% group_by summarise_each funs
#' @import NMF
NULL

#' Make heatmap illustrating the signel levels for clusters. For multiple samples, signal levels differ for the same cluster in different samples. Signal levels can be plotted easily using built-in R functions such as mean, median, and max.
#' 
#' @param data_agg Variable name of aggregated data.
#' @param signal_level Function used to find signal levels for visualization.
#' 
#' @return Saves the cluster signal level heatmap plot in the working directory.
#' @export

signalLevelHeatmap<-function(data_agg, signal_level = max){
  ClusterID<-NULL

  subpop_index<-c(min(data_agg[,1]):max(data_agg[,1]))
  sample_names<-unique(data_agg[,2])
  
  data_agg_sub<-data_agg[, -2]
  
  data_agg_intermediate<-data_agg_sub %>% group_by(ClusterID) %>% summarise_each(funs(signal_level))
  data_agg_intermediate<-data.frame(data_agg_intermediate)
  
  data_agg_new<- data_agg_intermediate[order(data_agg_intermediate$ClusterID),]

  M<-data_agg_new[,-c(1, ncol(data_agg_new))]
  #aheatmap(t(M))
  return(aheatmap(t(M)))
  
}
