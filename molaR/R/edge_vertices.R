#' Function for finding the edge vertices
#' 
#' Function will sort through all the vertices of the surface and find the ones which are on the edge. This will be needed for identifying which should be masked and not inlcuded in the calculation of the final DNE value.
#' @param plyFile a stanford PLY file
#'
#' 
#' edge_vertices()

edge_vertices <- function(plyFile) {
	Faces <- t(plyFile$it) # extract faces from ply object
	
	fnum <- length(Faces[,1])
	grump <- vector('list', fnum) ## new vector for holding ordered faces
	
	for (i in 1:fnum) {
		grump[[i]] <- unlist(sort(Faces[i,]))
	}
	
	edges <- vector('list', fnum) ## new vector for new created edge names
	
	for (i in 1:fnum) {
		grumps <- grump[[i]]
		
		G1 <- grumps[1]
		G2 <- grumps[2]
		G3 <- grumps[3]
		
		ED1 <- paste(G1, G2, sep='_')
		ED2 <- paste(G1, G3, sep='_')
		ED3 <- paste(G2, G3, sep='_')
		
		edges[[i]] <- paste(ED1, ED2, ED3, sep=',')
	}
	for (i in 1:fnum) {
		edges[[i]] <- unlist(strsplit(edges[[i]], ','))
	}
	
	verix <- matrix(0, fnum, 3) ## unpacking edge names
	for (i in 1:fnum) {
		gogo <- edges[[i]]
		verix[i,1] <- gogo[1]
		verix[i,2] <- gogo[2]
		verix[i,3] <- gogo[3]
	}
	
	string <- c(verix[,1], verix[,2], verix[,3]) ## repacking edge names
	
	edgeframe <- data.frame(names=string)
	
	grumpalump <- aggregate(edgeframe, list(edgeframe$names), FUN=length) ## counting number of repeated edges
	
	subgrump <- subset(grumpalump, grumpalump$names==1)
	
	uniques <- as.matrix(data.frame(uniques=subgrump$Group.1)) ## sorted out uniques
	
	onemore <- vector('list', length(uniques[,1]))
	for (i in 1:length(onemore)) {
		onemore[i] <- uniques[i,]
	}
	
	grumplestiltskin <- as.numeric(unlist(strsplit(unlist(onemore), '_'))) ## unpacking vertices
	
	tempframe <- data.frame(verts=grumplestiltskin)
	
	verts2 <- aggregate(tempframe, list(tempframe$verts), FUN=length)
	
	out <- as.vector(verts2$Group.1)
	
	return(out)
}
