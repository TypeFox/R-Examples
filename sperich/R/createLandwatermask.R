createLandwatermask <-
function(dataset.landwater,dimension, shift, resolution=1){
	#create grid
	landwatermask.nocoast <- matrix(0,dimension[1],dimension[2])

	#dataset available?
	if (!is.null(dataset.landwater)){
		#extract data
		long <- dataset.landwater$long
		lat <- dataset.landwater$lat
		landsum <- dataset.landwater$landsum

		#shift coordinates to origin
		long <- long - shift[1]
		lat <- lat - shift[2]

		#transformate coordinates to grid position
		long <- round(long / resolution) 
		lat <- round(lat / resolution)

		#insert data into grid
		for (i in 1:length(long)){
			if ((long[i]+1 > 0) && (long[i] < dimension[1]) && (lat[i]+1 > 0)&& (lat[i] < dimension[2]) && (landsum[i] == 0)){
				landwatermask.nocoast[long[i]+1,lat[i]+1] <- -1
			}
		}
	}

	return(landwatermask.nocoast)
}
