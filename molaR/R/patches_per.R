#' A function for patches within each face
#'
#' this gets some important information out of each patch
#' @param plyFile a stanford PLY file 
#' @param patch_details information on each patch
#' @param minimum_faces minimum number of faces in each counted patch
#' @param minimum_area minimum area for a patch to be counted
#' patches_per()


patches_per <- function(patch_details, plyFile, minimum_faces=3, minimum_area=0) {
	out <- list()
	
	patches_per <- matrix(0, nrow=length(names(patch_details)), ncol=1)
	rownames(patches_per) <- c(names(patch_details))
	if (minimum_area==0){
	for (i in 1:length(names(patch_details))) {
		A <- names(patch_details)[i]
		temp <- patch_details[[A]]
		 patches_per[A,1] <- length(which(temp[,1]>(minimum_faces-1)))
	}} else {
		MinAreaPercentage <- sum(plyFile$Face_Areas)*minimum_area
		for (i in 1:length(names(patch_details))) {
		A <- names(patch_details)[i]
		temp <- patch_details[[A]]
		 patches_per[A,1] <- length(which(temp[,2]>MinAreaPercentage))
	}
	}
	out[['total patches']] <- sum(patches_per[,1])
	out[['directions']] <- patches_per
	return(out)
}