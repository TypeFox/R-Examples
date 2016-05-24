# function to summarize the PAFitData object  2015-3-11 Thong Pham
summary.PAFitData <- function(object,...){
  print(object)
  cat("Type of network: ",object$net_type[1],"\n");
  cat("Number of nodes in the final network: ",length(object$final_deg),"\n")
  cat("Number of edges in the final network: ",sum(object$final_deg),"\n")
  cat("Number of new nodes: ",length(object$final_deg) - object$initial_nodes,"\n")
  cat("Number of new edges: ",sum(object$Sum_m_k), "\n")
  cat("Numbef of time-steps:",  object$T,"\n");
  cat("Maximum degree: ",object$deg.max,"\n");
  cat("Number of bins:",object$G,"\n");
  cat("Threshold of the number of new edges acquired: ", object$deg_thresh,"\n")
  cat("Number of nodes satisfied the threshold: ",length(object$f_position),"\n")
}