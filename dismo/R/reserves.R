# Author: Robert J. Hijmans, r.hijmans@gmail.com
# 2009
# Version 0.1
# Licence GPL3


.rebelo <- function(xy) {
	xy[is.na(xy)] <- 0
	xy <- round(xy)
	xy[xy<0] <- 0
	xy[xy>0] <- 1
	nspecies <- nrow(xy)
	nsites <- ncol(xy)
	res <- matrix(ncol=2, nrow=nsites)
	for (i in 1:nsites) {
		sitesppcount <- colSums(xy)
		nsp <- max(sitesppcount)
		if (nsp == 0) {break}
		selsite <- which(sitesppcount == nsp)[1]
		res[i,1] <- selsite
		res[i,2] <- nsp
		delspp <- as.vector(which(xy[,selsite]==1))
		xy[delspp,] <- 0
	}
	colnames(res) <- c("Site", "nSpecies")
	return(res[1:(i-1),])
}

	