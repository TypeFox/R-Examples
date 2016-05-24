colorVectorIsNull <-
function(data_matrix){
	design_matrix <- matrix(1,nrow(data_matrix),1)
	rownames(design_matrix) <- rownames(data_matrix)
	color_vector <- createColorVectorsByDesign(design_matrix)
	return(color_vector)
}
