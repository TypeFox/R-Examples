nucDiag <- function(DNAbin, sppVector){
	DNAbin <- as.matrix(DNAbin)
	inform <- seg.sites(DNAbin)
	sppSeqs <- lapply(unique(sppVector), function(x) which(sppVector == x))
	
	siteCheck <- function(spp, site){
		res <- as.character(DNAbin[spp, site]) %in% as.character(DNAbin[-spp, site])
		#A 'res' of TRUE means that the nucleotide in the sp. is also present in the rest of the spp.
		res <- as.logical(sum(as.numeric(res)))
		res
	}
	li <- list()
	for(i in 1:length(sppSeqs)){
		li[[i]] <- NA
		for(j in 1:length(inform)){
			li[[i]][j] <- siteCheck(sppSeqs[[i]], inform[j])
		}
	}
	out <- lapply(li, function(x) inform[which(!x)])
	names(out) <- unique(sppVector)
	out
}
