species.richness.cv <- 
function(dataset.all.species, landwatermask, fold=5, loocv.limit=10, 
		distances=3:10, weight=0.5, dimension, shift, resolution=1, 
		upperbound, all.species=-1, silent=TRUE, do.parallel=FALSE){
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
	species.richness.weighted.cv <- matrix(0,dimension[1],dimension[2])

	#iterate parallel about all species
	requireNamespace("foreach")
	if (do.parallel){
		if(!foreach::getDoParRegistered()) {
			if (!silent) {cat("No parallel backend detected! Problem will be solved sequential.\n",sep="")}
			foreach::registerDoSEQ()
		} else {
			if (!silent) {cat("Parallel backend detected.\n",sep="")}
		}

		species.richness.weighted.cv <- foreach::foreach(species=all.species, .combine="+", .inorder=FALSE) %dopar% {
			#create grid
			species.richness.weighted.one.species <- matrix(0,dimension[1],dimension[2])
			
			#initialize temporary array
			species.range.all.subs <- array(0, dim=c(length(distances),dimension[1],dimension[2]))

			#extract datasets of one species out of database
			dataset.one.species <- extract.species(dataset.all.species, species)
			number.of.occurrences <- dim(dataset.one.species)[1]

			#more than two occurrences needed
			if (number.of.occurrences > 2){
				subsamples <- generate.subsamples(number.of.occurrences, fold, loocv.limit)

				#iterate about distances
				for (distance in distances){
					#iterate about all subsamples
					for (subsample.id in 1:dim(subsamples)[1]){
						subsample <- subsamples[subsample.id,]
						subsample <- subsample[which(subsample != 0)]
						dataset.one.subsample <- dataset.one.species[subsample,]

						#calculate species range	
						species.range.sub <- species.range(dataset.one.subsample, distance, dimension, 
									shift, resolution, landwatermask, upperbound, cross.validation=TRUE)
					
						#sum over all subsamples
						species.range.all.subs[which(distance == distances),,] <- species.range.all.subs[which(distance == distances),,] + species.range.sub
					}				

					#divide through number of subsamples
					species.range.sub.tmp <- species.range.all.subs[which(distance == distances),,] / matrix(dim(subsamples)[1],dimension[1],dimension[2])
					species.range.sub.tmp[which(is.na(species.range.sub.tmp)==TRUE)] <- 0
					species.range.all.subs[which(distance == distances),,] <- species.range.sub.tmp				
	
					if (which(distance==distances)==1){
						species.richness.weighted.one.species <- species.range.all.subs[1,,]
					} else {
						species.richness.weighted.one.species <- species.richness.weighted.one.species + 
							(distance^(-weight) * (species.range.all.subs[which(distance == distances),,] - species.range.all.subs[which(distance == distances)-1,,]))
					}
	
				}
				
				# output message
				if (!silent){
					cat(rep("\b", nchar(message)),sep="")
					message <- paste("Species ",which(species==all.species)," of ",number.of.species," done!", sep="")
					cat(message)
					flush.console()
				}
			}
			# foreach return statement
			return(species.richness.weighted.one.species)	
		}

	} else {
	#iterate about all species
		for (species in all.species){
			#create grid
			species.richness.weighted.one.species <- matrix(0,dimension[1],dimension[2])

			#initialize temporary array
			species.range.all.subs <- array(0, dim=c(length(distances),dimension[1],dimension[2]))

			#extract datasets of one species out of database
			dataset.one.species <- extract.species(dataset.all.species, species)

			number.of.occurrences <- dim(dataset.one.species)[1]

			#more than two occurrences needed
			if (number.of.occurrences > 2){
				subsamples <- generate.subsamples(number.of.occurrences, fold, loocv.limit)

				#iterate about distances
				for (distance in distances){
					#iterate about all subsamples
					for (subsample.id in 1:dim(subsamples)[1]){
						subsample <- subsamples[subsample.id,]
						subsample <- subsample[which(subsample != 0)]
						dataset.one.subsample <- dataset.one.species[subsample,]

						#calculate species range	
						species.range.sub <- species.range(dataset.one.subsample, distance, dimension, 
									shift, resolution, landwatermask, upperbound, cross.validation=TRUE)
					
						#sum over all subsamples
						species.range.all.subs[which(distance == distances),,] <- species.range.all.subs[which(distance == distances),,] + species.range.sub
					}				

					#divide through number of subsamples
					species.range.sub.tmp <- species.range.all.subs[which(distance == distances),,] / matrix(dim(subsamples)[1],dimension[1],dimension[2])
					species.range.sub.tmp[which(is.na(species.range.sub.tmp)==TRUE)] <- 0
					species.range.all.subs[which(distance == distances),,] <- species.range.sub.tmp				
	
					if (which(distance==distances)==1){
						species.richness.weighted.one.species <- species.range.all.subs[1,,]
					} else {
						species.richness.weighted.one.species <- species.richness.weighted.one.species + 
							(distance^(-weight) * (species.range.all.subs[which(distance == distances),,] - species.range.all.subs[which(distance == distances)-1,,]))
					}
	
				}

				#sum over all species
				species.richness.weighted.cv <- species.richness.weighted.cv + species.richness.weighted.one.species
			
				if (!silent){
					cat(rep("\b", nchar(message)),sep="")
					message <- paste("Species ",which(species==all.species)," of ",number.of.species," done!", sep="")
					cat(message)
					flush.console()
				}
			}
		}		
	}

	if (!silent)
		cat("\n")

	return(species.richness.weighted.cv)
}
