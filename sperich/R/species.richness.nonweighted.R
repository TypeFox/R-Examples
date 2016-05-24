species.richness.nonweighted <- 
function(dataset.all.species, landwatermask, distance=10, 
		dimension, shift, resolution=1, upperbound, 
		all.species=-1, silent=TRUE, do.parallel=FALSE){
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
	message <- ""

	#create grid
	species.richness.noweight <- matrix(0, dimension[1], dimension[2])

	requireNamespace("foreach")
	if (do.parallel){
		if(!foreach::getDoParRegistered()) {
			if (!silent) {cat("No parallel backend detected! Problem will be solved sequential.\n",sep="")}
			foreach::registerDoSEQ()
		} else {
			if (!silent) {cat("Parallel backend detected.\n",sep="")}
		}
	
		species.richness.noweight <- foreach::foreach(species=all.species, .combine="+", .inorder=FALSE) %dopar% {
			dataset.one.species <- extract.species(dataset.all.species, species)
			species.range.d <- species.range(dataset.one.species, distance, 
						dimension, shift, resolution, landwatermask, 
						upperbound)
			if (!silent){
				cat(rep("\b", nchar(message)),sep="")
				message <- paste("Species ",which(species==all.species)," of ",number.of.species," done!", sep="")
				cat(message)
				flush.console()
			}	
			return(species.range.d)
		}
	} else {
		for (species in all.species){
			dataset.one.species <- extract.species(dataset.all.species, species)
			species.range.d <- species.range(dataset.one.species, distance, 
						dimension, shift, resolution, landwatermask, 
						upperbound)
			#sum over all species
			species.richness.noweight <- species.richness.noweight + species.range.d
		
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

	return(species.richness.noweight)
}
