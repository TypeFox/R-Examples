
queryGISD <- function(species) {
	
	species <- gsub(' ', '_', species)
	
	if (species %in% names(invasiveDB[[1]])) {	
		native <- invasiveDB[['nativeRanges']][[species]]
		alien <- invasiveDB[['alienRanges']][[species]]
	} else {
		native <- 'not found'
		alien <- 'not found'
	}
	return(list(species=species, native=native, alien=alien))	
}

