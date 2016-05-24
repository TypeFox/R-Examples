#' Function for gathering the patches for each direction
#'
#' This function will gather the patches in each of the 8 bins and ready it for
#' patches_for_each_direction()
#' @param indexed_pairs Pairs of touching faces

patches_for_each_direction <- function(indexed_pairs){
clusterlist <- list()

for (i in names(indexed_pairs)){
	direction <- indexed_pairs[[as.character(i)]]
	clusterlist[[as.character(i)]] <- clustered_patches(direction)
	
}

return(clusterlist)
}
