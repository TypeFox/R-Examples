sppDist <-
function(distobj, sppVector){
	distobj <- as.matrix(distobj)
	attr(distobj, "dimnames")[[1]] <- sppVector
	taxa <- unique(sppVector)
	intra <- list()
	inter <- list()
	for(i in 1:length(taxa)){
		for(j in 1:length(taxa)){
			sppMat <- distobj[which(dimnames(distobj)[[1]] == taxa[i]), which(dimnames(distobj)[[1]] == taxa[j])]
			if(taxa[i] == taxa[j]) intra[[length(intra)+1]] <- sppMat[lower.tri(sppMat)] else inter[[length(inter)+1]] <- sppMat
			}
	}
list(inter = unlist(inter), intra = unlist(intra))
}

