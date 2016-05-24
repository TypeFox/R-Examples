cophen <- function(phy){
	n=Ntip(phy)
	out <- .Call("cophen", tree=list(
								   ROOT = as.integer(n+1),
								   MAXNODE = as.integer(max(phy$edge[,1])),
								   ENDOFCLADE = as.integer(dim(phy$edge)[1]),
								   ANC = as.integer(phy$edge[,1]),
								   DES = as.integer(phy$edge[,2]),
								   EDGES = as.double(c(phy$edge.length,0)),
								   COPHEN = as.double(array(matrix(0, n, n)))),
				 PACKAGE = "spacodiR")
	cc=matrix(out$COPHEN,nrow=n,byrow=FALSE)
	rownames(cc)<-colnames(cc)<-phy$tip.label
	return(cc)
}