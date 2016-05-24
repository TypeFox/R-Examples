edgeNotValid <- 
function(grid, point.a, point.b, landwatermask, upperbound){
	#save old grid
	old_grid <- grid
	#add edge to grid
	new_grid <- add.Edges(grid, point.a,point.b)
	#extract changes
	edges <- new_grid - old_grid
	#extract position of changes
	edges.position <- which(edges > 0)
	
	#check if new edge crosses water or high altitude
	if ((length(which(landwatermask[edges.position]==-1)) > 0)||(length(which(landwatermask[edges.position] >= upperbound)) > 0)){
		return(TRUE)
	} else {
		return(FALSE)
	}
}