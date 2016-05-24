fill.Areas <- 
function(grid, landwatermask){
	#speed up iteration (reducing points)
	maxlong <- dim(grid)[1]
	maxlat <- dim(grid)[2]
	minlong <- minlat <- 0
		
	##reduce x-steps
	for (i in 1:dim(grid)[1]){
		if (length(which(grid[i,] > 0)) > 0){
			maxlong <- i
		}
		if ((length(which(grid[i,] > 0)) > 0)&&(minlong==0)){
			minlong <- i
		}
	}
	if (minlong == maxlong){
		return(grid)
	}

	##reduce y-steps
	for (i in 1:dim(grid)[2]){
		if (length(which(grid[,i] > 0)) > 0){
			maxlat <- i
		}
		if ((length(which(grid[,i] > 0)) > 0)&&(minlat==0)){
			minlat <- i
		}
	}
	if (minlat == maxlat){
		return(grid)
	}
	
	#floodfill
	flaggrid <- newgrid <- grid[minlong:maxlong, minlat:maxlat]
	landwatermask.newdim <- landwatermask[minlong:maxlong, minlat:maxlat]
	queue <- list()

	for (x in 1:(dim(newgrid)[1])){
		queue[[length(queue)+1]] <- c(x,1)
		queue[[length(queue)+1]] <- c(x,dim(newgrid)[2])
	}
	
	for (y in 1:(dim(newgrid)[2])){
		queue[[length(queue)+1]] <- c(1,y)
		queue[[length(queue)+1]] <- c(dim(newgrid)[1],y)

	}

	while(length(queue) > 0){
		point <- queue[[1]]
		queue[[1]] <- NULL
		
		if (flaggrid[point[1],point[2]] <= 0){
			flaggrid[point[1],point[2]] <- 1

			if ((point[2]+1) <= dim(newgrid)[2]){
				queue[[length(queue)+1]] <- c(point[1],point[2]+1)
			}
			if ((point[2]-1) > 0){
				queue[[length(queue)+1]] <- c(point[1],point[2]-1)
			}
			if ((point[1]+1) <= dim(newgrid)[1]){
				queue[[length(queue)+1]] <- c(point[1]+1,point[2])
			}
			if ((point[1]-1) > 0){
				queue[[length(queue)+1]] <- c(point[1]-1,point[2])
			}
		}
	}

	check <- which(flaggrid < 1)	
	newgrid[check[which(landwatermask.newdim[check] >= 0)]] <- 1

	grid[minlong:maxlong, minlat:maxlat] <- newgrid

	return(grid)
}
