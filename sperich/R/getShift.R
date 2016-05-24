getShift <- 
function(dataset.all.species){
	return(c(min(dataset.all.species$long), min(dataset.all.species$lat)))
}
