createColorVectorsByDesign <-
function(design_matrix,hsv=TRUE,offset=NULL){

	if(hsv){
		group_colors <- prettyGraphsHSVColorSelection(n.colors=ncol(design_matrix),offset=offset)
	}else{
		group_colors <- prettyGraphsColorSelection(n.colors=ncol(design_matrix),offset=offset)	
	}
	rownames(group_colors)<-colnames(design_matrix)

	arr.ind <- which(design_matrix==1,arr.ind=TRUE)
	arr.ind <- arr.ind[order(arr.ind[,1]),]
	observation_colors <- as.matrix(group_colors[arr.ind[,2],1])
	rownames(observation_colors) <- rownames(design_matrix)
	return(list(oc=observation_colors,gc=group_colors))
}
