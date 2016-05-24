getDimension <- 
function(dataset.all.species, resolution=1){
	#extract coordinates
	long <- dataset.all.species$long
	lat <- dataset.all.species$lat
	
	#shift coordinates to origin
	long <- long - min(long)
	lat <- lat - min(lat)

	#transformate coordinates to grid position
	long <- round(long / resolution )
	lat <- round(lat / resolution )
	
	return(c(max(long) + 1,max(lat)+1))
}
