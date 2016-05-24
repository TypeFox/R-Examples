predgrid <-
function(dataXY, map=NULL, nrow=100, ncol=100) {
	XYcol = 1:2 + (length(dataXY) > 2)		# check for outcome column
	dataXrange=range(dataXY[,XYcol[1]])
	dataYrange=range(dataXY[,XYcol[2]])
	# Creates a rectangular grid using range of data X and Y
	grid=as.data.frame(expand.grid(X=seq(dataXrange[1],dataXrange[2],len=ncol),
								   Y=seq(dataYrange[1],dataYrange[2],len=nrow)))
	# If input data geolocation columns have names assign them to the output
	if (!is.null(names(dataXY))) names(grid) = names(dataXY)[XYcol]
	# If map is provided, grid is clipped using trimdata function
	if (!is.null(map)) grid = trimdata(grid,map) 					
	return(grid)	
}
