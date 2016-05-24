species.richness <- 
function(dataset.all.species, landwatermask, distances=2:10, weight=0.5, 
		dimension, shift, resolution=1, upperbound,  
		all.species=-1, silent=TRUE, do.parallel=FALSE){
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
	
	#adjust distances to include baseline distance
	distances <- c(min(distances)-1, distances)

	message <- ""

	#create grids
	species.richness.weighted <- matrix(0, dimension[1], dimension[2])
	species.range.distance <- array(0, dim=c(length(distances),dimension[1], dimension[2]))

	requireNamespace("foreach")
	if (do.parallel){
		# check if parallel backend is available
		if(!foreach::getDoParRegistered()) {
			if (!silent) {cat("No parallel backend detected! Problem will be solved sequential.\n",sep="")}
			foreach::registerDoSEQ()
		} else {
			if (!silent) {cat("Parallel backend detected.\n",sep="")}
		}

		# iterate parallel or sequential with foreach-loop
		species.richness.weighted <- foreach::foreach(species=all.species, .combine="+", .inorder=FALSE) %dopar% {
			# create species grid
			species.richness.weighted.one.species <- matrix(0, dimension[1], dimension[2])
			# extract species
			dataset.one.species <- extract.species(dataset.all.species, species)
			# iterate about distances
			for (distance in distances){
				species.range.distance[which(distance == distances),,] <- species.range(dataset.one.species, distance, 
											dimension, shift, resolution, landwatermask, 
											upperbound)
	
				if (which(distance==distances)==1){
					species.richness.weighted.one.species <- species.range.distance[1,,]
				} else {
					species.richness.weighted.one.species <- species.richness.weighted.one.species + 
					(distance^(-weight) * (species.range.distance[which(distance == distances),,] - species.range.distance[which(distance == distances)-1,,]))
				}
			}
			# output message
			if (!silent){
				cat(rep("\b", nchar(message)),sep="")
				message <- paste("Species ",which(species==all.species)," of ",number.of.species," done!", sep="")
				cat(message)
				flush.console()
			}
			# foreach return statement
			return(species.richness.weighted.one.species)
		}
	} else {
		# iterate sequentially (with for-loop)
		for (species in all.species){
			# create species grid
			species.richness.weighted.one.species <- matrix(0, dimension[1], dimension[2])
			# extract species
			dataset.one.species <- extract.species(dataset.all.species, species)
			# iterate about distances
			for (distance in distances){
				species.range.distance[which(distance == distances),,] <- species.range(dataset.one.species, distance, 
											dimension, shift, resolution, landwatermask, 
											upperbound)
	
				if (which(distance==distances)==1){
					species.richness.weighted.one.species <- species.range.distance[1,,]
				} else {
					species.richness.weighted.one.species <- species.richness.weighted.one.species + 
					(distance^(-weight) * (species.range.distance[which(distance == distances),,] - species.range.distance[which(distance == distances)-1,,]))
				}
	
			}
			
			#sum over all species
			species.richness.weighted <- species.richness.weighted + species.richness.weighted.one.species
			
			# output message
			if (!silent){
				cat(rep("\b", nchar(message)),sep="")
				message <- paste("Species ",which(species==all.species)," of ",number.of.species," done!", sep="")
				cat(message)
				flush.console()
			}	
		}
	}
	
	if (!silent){
		cat("\n")
	}

	return(species.richness.weighted)
}
