species.range <- 
function(dataset.one.species, distance, dimension, shift, 
		resolution=1, landwatermask, upperbound=5, 
		cross.validation=FALSE){	
	#create grid
	grid <- matrix(0,dimension[1],dimension[2])

	#add points
	grid <- data.into.Grid(dataset.one.species, dimension, shift, resolution)

	#points into list
	points <- which(grid > 0)
	points.xy <- list()
	
	for (i in 1:length(points)){
		points.xy[[i]] <- c(ifelse((points[i] %% dimension[1]) == 0, dimension[1], points[i] %% dimension[1]), 
					ceiling(points[i]/dimension[1]))
	}

	if (cross.validation){
		#test for neighbours
		points.valid <- list()
		points.xy.old <- points.xy

		while(length(points.valid) != length(points.xy.old)){
			points.valid <- list()
			for (i in 1:length(points.xy)){
				point <- points.xy[[i]]
				neighbour.found <- 0
				for (j in 1:length(points.xy)){
					if (i != j){
						if (getDistance(point,points.xy[[j]], resolution) <= distance){
							neighbour.found <- neighbour.found + 1
						}
					}
					if (neighbour.found >= 2){
						break
					}
				}	
				if (neighbour.found >= 2){
					points.valid[[length(points.valid)+1]] <- point
				}
	 		}

			if (length(points.valid) < 3){
				return(matrix(0,dimension[1],dimension[2]))
			}	

			points.xy.old <- points.xy
			points.xy <- points.valid
		}

		#create new grid without invalid points
		grid <- matrix(0,dimension[1],dimension[2])
		for (point in points.xy){
			grid[point[1],point[2]] <- 1
		}
	}

	#add edges
	for (i in 1:length(points.xy)){
		for (j in 1:length(points.xy)){
			if ((i!=j)&&(getDistance(points.xy[[i]],points.xy[[j]], resolution) <= distance )&&(!edgeNotValid(grid, points.xy[[i]],points.xy[[j]], landwatermask, upperbound))){
				grid <- add.Edges(grid, points.xy[[i]],points.xy[[j]])
			}
		}
	}

	#add areas
	grid <- fill.Areas(grid, landwatermask)

	return(grid)
}
