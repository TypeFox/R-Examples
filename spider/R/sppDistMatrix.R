sppDistMatrix <-
function(distobj, sppVector){
	dat <- as.matrix(distobj)
	attr(dat, "dimnames")[[1]] <- sppVector
	taxa <- unique(sppVector)
	pair.mat <- matrix(data = NA, nrow = length(taxa), ncol = length(taxa), dimnames = list(one = taxa,two = taxa))
	for(i in 1:length(taxa)){
		for(j in 1:length(taxa)){
			sppMat <- dat[which(dimnames(dat)[[1]] == taxa[i]), which(dimnames(dat)[[1]] == taxa[j])]
			if(taxa[i] == taxa[j]) pair.mat[i,j] <- mean(sppMat[lower.tri(sppMat)]) else pair.mat[i,j] <- mean(sppMat)
		}
	}
pair.mat
}

