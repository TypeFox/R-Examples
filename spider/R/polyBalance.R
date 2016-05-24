polyBalance <- function(phy){
	nd <- node.depth(phy)
	tab <- table(phy$edge[,1])
	ind <- as.numeric(names(tab[which(tab == 2)]))
	tips <- nd[phy$edge[phy$edge[, 1] %in% ind, 2]]
	nodes <- phy$edge[phy$edge[, 1] %in% ind, 1]
	mat <- cbind(nodes, tips)
	NNodes <- unique(phy$edge[,1])
	matInd <- unique(mat[,1])
	outMat <- matrix(NA, ncol = 2, nrow = length(NNodes))
	dimnames(outMat)[[1]] <- NNodes
	for(i in 1:length(matInd)){
		outMat[NNodes == matInd[i],] <- mat[mat[,1] == matInd[i],2]
	}
	
	outMat
}
