extract.species <- 
function(dataset.all.species, species.number){
	select <- which(dataset.all.species$speciesID==species.number)
	dataset.one.species <- data.frame(long=dataset.all.species$long[select], lat=dataset.all.species$lat[select], 
						speciesID=dataset.all.species$speciesID[select])
	return(dataset.one.species)
}