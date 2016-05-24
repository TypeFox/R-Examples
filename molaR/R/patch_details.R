#' Function for gathering patch details for each Orientation patch
#' 
#' This function does some simple math to lets us know about the patches
#' @param plyFile a stanford PLY file 
#' @param clusterlist a list of faces in the cluster
#' patch_details()


patch_details <- function(clusterlist, plyFile) {
cluster_list_details <- list()

for (i in names(clusterlist)) {
	
	temp <- clusterlist[[as.character(i)]]
	
	cluster_stats <- matrix(0, nrow=length(names(temp)), ncol=2)
	colnames(cluster_stats) <- c('Faces', 'Total Area')
	rownames(cluster_stats) <- c(paste(as.character(i), names(temp), sep=":"))
	
	
	for (j in 1:length(names(temp))) {
		cluster_stats[j,1] <- length(temp[[j]])
		cluster_stats[j,2] <- sum(plyFile$Face_Areas[as.numeric(unlist(temp[[j]]))])
	}
	cluster_list_details[[as.character(i)]] <- cluster_stats
}
return(cluster_list_details)

}
