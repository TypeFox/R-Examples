data.into.Grid <-
function(dataset.one.species, dimension, shift, resolution=1){
	#create grid
	grid <- matrix(0,dimension[1],dimension[2])

	#extract coordinates
	long <- dataset.one.species$long
	lat <- dataset.one.species$lat


	#shift coordinates to origin
	long <- long - shift[1]
	lat <- lat - shift[2]

	#transformate coordinates to grid position
	long <- round(long / resolution) 
	lat <- round(lat / resolution)

	#insert data into grid
	for (i in 1:length(long)){
		grid[long[i]+1,lat[i]+1] <- 1	
	}
	return(grid)
}
