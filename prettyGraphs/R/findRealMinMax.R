findRealMinMax <-
function(mat1,mat2,axis1,axis2){
	track_min_max_matrix <- rbind(mat1,mat2)
	min_max_list <- list(minx=min(c(track_min_max_matrix[,axis1],0)),miny=min(c(track_min_max_matrix[,axis2],0)),maxx=max(c(track_min_max_matrix[,axis1],0)),maxy=max(c(track_min_max_matrix[,axis2],0)))
	return(min_max_list)
}
