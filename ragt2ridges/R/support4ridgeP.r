support4ridgeP <- function(adjMat=NULL, nNodes=NULL, zeros=NULL, verbose=FALSE){
	#####################################################################################################
	# 
	# DESCRIPTION:
	# The zero pattern of an adjancency matrix is converted into its junction tree (as a list of its cliques and separators). 
	# When a graph with nondecomposable support is provided, it is triangulated.
	#
	# ARGUMENTS:
	# -> adjMat     : Adjacency matrix of an undirected graph.
	# -> nNodes     : Number of nodes of the graph.
	# -> zeros      : Matrix with indices of entries of the adjacency matrix that are constrained to zero. The matrix 
	#                 comprises two columns, each row corresponding to an entry of the adjacency matrix. The first 
	#                 column contains the row indices and the second the column indices. The specified graph
	#                 should be undirected and decomposable. If not, it is symmetrized and triangulated.
	# -> verbose	: Logical indicating whether progress should be reported.
	# 
	#
	# DEPENDENCIES: 
	# require("igraph")          # functions from package : graph.adjancency, igraph.to.graphNEL
	# require("gRbase")          # functions from package : triangulate
	# require("RGBL")            # functions from package : is.triangulated
	# require("graph")           # functions from package : 'graphAM'-class 
	#
	#####################################################################################################

	# iput checks
	if (is.null(adjMat) & is.null(nNodes) & is.null(zeros)){ stop("Support not sufficiently specified.") }
	if (is.null(adjMat) & is.null(nNodes)){ stop("Support not sufficiently specified.") }
	if (is.null(adjMat) & is.null(zeros)){ stop("Support not sufficiently specified.") }
	if (!is.null(adjMat)){
		if (!is.matrix(adjMat)){ stop("adjMat should be a matrix.") }
		if (nrow(adjMat) != ncol(adjMat)){ stop("adjMat should be square matrix.") }
		if (!all(sort(unique(as.numeric(adjMat))) == c(0,1))){ stop("Elements of adjMat ill-specified.") }
		if (!isSymmetric(adjMat)){ stop("adjMat does not correspond to an undirect graph.") }
		if (!all(diag(adjMat) == 0)){ stop("adjMat contains nonzero diagonal element.") }
	}
	if (!is.null(nNodes)){
		if (class(nNodes) != "numeric" & class(nNodes) != "integer"){ stop("Input (nNodes) is of wrong class") }
		if (length(nNodes) != 1){ stop("Input (nNodes) must be a scalar") }
	      	if (!.is.int(nNodes)){ stop("Input (nNodes) should be an integer") }
		if (nNodes <= 0){ stop("Input (nNodes) must be strictly positive") }
	}
	if (!is.null(zeros)){
		if (!is.null(zeros) & as.character(class(zeros)) != "matrix"){ stop("Input (zeros) is of wrong class.") }    
		if (!is.null(zeros)){ if(ncol(zeros) != 2){ stop("Wrong dimensions of the (zeros).") } } 
		if (!is.null(zeros)){ zeros <- zeros[order(zeros[,2], zeros[,1]),] }
	}
	if (as.character(class(verbose)) != "logical"){ stop("Input (verbose) is of wrong class.") }

	# convert zero pattern into adjacency matrix
	if (!is.null(zeros) & !is.null(nNodes) & is.null(adjMat)){
		adjMat <- matrix(1, nNodes, nNodes)
		adjMat[zeros] <- 0
		adjMat[cbind(zeros[,2], zeros[,1])] <- 0
		diag(adjMat) <- 0
	}

	# convert adjacency into graphNel object
	G <- igraph.to.graphNEL(graph.adjacency(adjMat, mode="undirected"))

	# is graph complete?
	is.complete(G)

	# check whether a chordal support has been provided
	addedEdges <- matrix(nrow=0, ncol=2)
	if (!is.triangulated(G)){ 
		if(verbose){ cat("-> provided zero pattern not chordal   : support is triangulated,", "\n") }
		nEdgesOld <- numEdges(G)
		G <- triangulate(G)
		nEdgesNew <- numEdges(G)
		addedEdges <- which(as(G, "graphAM")@adjMat - adjMat == 1, arr.ind=TRUE)
		adjMat <- as(G, "graphAM")@adjMat
		diag(adjMat) <- 1
		zeros <- which(adjMat == 0, arr.ind=TRUE)
		if(verbose){ cat(paste("->                                       ", nEdgesNew - nEdgesOld, " extra nonzeros", sep=""), "\n") }
	}
	if (is.null(zeros)){
		diag(adjMat) <- 1
		zeros <- which(adjMat == 0, arr.ind=TRUE)
	}

	# decompose in cliques and separators
	decomp <- rip(G)
	Cliques <- decomp$cliques
	Separators <- decomp$separators

	# convert index vectors from character to integer
	for (k in 1:length(Cliques)){ Cliques[[k]] <- as.integer(Cliques[[k]]) }
	for (k in 1:length(Separators)){ Separators[[k]] <- as.integer(Separators[[k]]) }

	return(list(zeros=zeros, cliques=Cliques, separators=Separators, addedEdges=addedEdges))
}
