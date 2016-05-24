getNarrowEndemics <- function(dataset.all.species, all.species, 
			narrow.endemic.limit, dimension, shift, 
			resolution){
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

	# function to check narrow endemism of single species
	isNarrowEndemic <- function(dataset.one.species, dimension, shift, resolution){
		# check number of occurrences
		if (dim(dataset.one.species)[1] <= narrow.endemic.limit){
			#check distance between occurrences
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
			
			for (i in 1:length(points.xy)){
				point <- points.xy[[i]]
				for (j in 1:length(points.xy)){
					if (i != j){
						if (getDistance(point,points.xy[[j]], resolution) > narrow.endemic.limit){
							return(FALSE)
						}
					}
				}
			}
			return(TRUE)

		}
		return(FALSE)
	}

	all.narrow.endemics <- c()

	# check all species
	for (species in all.species){
		dataset.one.species <- extract.species(dataset.all.species, species)
		if (isNarrowEndemic(dataset.one.species, dimension, shift, resolution)){
			all.narrow.endemics <- c(all.narrow.endemics, species)
		}
	}

	return(all.narrow.endemics)
}