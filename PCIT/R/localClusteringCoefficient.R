localClusteringCoefficient <- function(adj) {
	result <- .Fortran( "localCC", adj=as.single(adj), nGenes=as.integer(nrow(adj)), cc=as.single(rep(0, length=nrow(adj))), E=as.single(rep(0,nrow(adj))), k=as.integer(rep(0,nrow(adj))) )
	
	return(as.vector(result$cc))
}
