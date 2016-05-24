createHeightmask <-
function(dataset.height,dimension, shift, resolution=1){
	#create grid
	height.matrix <- matrix(0,dimension[1],dimension[2])

	#dataset available?
	if (!is.null(dataset.height)){
		#extract data
		long <- dataset.height$long
		lat <- dataset.height$lat
		height <- dataset.height$height

		#shift data to origin
		long <- long - shift[1]
		lat <- lat - shift[2]

		#transformate coordinates to grid position
		long <- round(long / resolution) 
		lat <- round(lat / resolution)

		#insert data into grid
		for (i in 1:length(long)){
			if (((long[i]+1) > 0) && (long[i] < dimension[1]) && ((lat[i]+1) > 0)&& (lat[i] < dimension[2])){
				if (height.matrix[long[i]+1,lat[i]+1] < height[i]){
					height.matrix[long[i]+1,lat[i]+1] <- height[i]
				}
			}
		}
	}

	return(height.matrix)
}
