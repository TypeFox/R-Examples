determineAxesPosition <-
function(min_max_list){
#	if(is.null(supplementary_matrix)){
#		minMaxList <- minmaxHelper(data_matrix,data_matrix,x_axis,y_axis)
#	}else{
#		minMaxList <- minmaxHelper(data_matrix,supplementary_matrix,x_axis,y_axis)
#	}
	minx <- min_max_list$minx
	miny <- min_max_list$miny
	maxx <- min_max_list$maxx
	maxy <- min_max_list$maxy
	x_axis_y_pos <- c(0,0)
	y_axis_x_pos <- c(0,0)
	y_axis_y_pos <- c(miny,maxy)
	x_axis_x_pos <- c(minx,maxx)
	return(list(xx=x_axis_x_pos,xy=x_axis_y_pos,yy=y_axis_y_pos,yx=y_axis_x_pos))
}
