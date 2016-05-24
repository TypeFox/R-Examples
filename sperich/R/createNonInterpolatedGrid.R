createNonInterpolatedGrid <-
function(dataset.all.species, dimension, shift, resolution=1, all.species=-1){
	#check species
	if (all.species[1]==-1){
		all.species <- unique(dataset.all.species$speciesID)
	} else {
		all.species.tmp <- c()
		for (species in all.species){
			if (length(which(dataset.all.species$speciesID == species)==TRUE) > 0){
				all.species.tmp <- c(all.species.tmp, species)
			}
		}
		all.species <- all.species.tmp
	}
	number.of.species <- length(all.species)
	
	#create grid
	grid <- matrix(0,dimension[1],dimension[2])

	for (species in all.species){
		dataset.one.species <- extract.species(dataset.all.species, species)
		grid <- grid + data.into.Grid(dataset.one.species, dimension, shift, resolution)
	}

	return(grid)
}