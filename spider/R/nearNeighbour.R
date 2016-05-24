nearNeighbour <- function(distobj, sppVector, names = FALSE){
	distobj <- as.matrix(distobj)
	diag(distobj) <- NA
	aa <- apply(distobj, MARGIN=2, FUN=function(x) which(x == min(x, na.rm = TRUE)))
	bb <- sapply(aa, function(x) names(sort(table(sppVector[x]), decreasing=TRUE)[1]))
	if(names) as.vector(bb) else as.vector(bb == sppVector)
}
