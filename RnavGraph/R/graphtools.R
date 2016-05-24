setMethod(
		f = "linegraph",
		signature = "graph",
		definition = function(graph,sep=":"){
			
			if(isDirected(graph)){
				stop("Graph must not be directed an graph is directed.")
			}else{
				
				## TODO test with duplicates=FALSE
				## from-to matrix
				eM <- edgeMatrix(as(graph,"graphNEL"), duplicates=FALSE)
				
				## node names of new graph
				V <- paste(nodes(graph)[eM[1,]],nodes(graph)[eM[2,]],sep=sep)
				
				new_eM <- NULL
				
				n <- dim(eM)[2]
				for(i in 1:n){
					for(ii in i:n){
						if(i!=ii){
							if(any(match(eM[,i],eM[,ii],nomatch=FALSE))){
								new_eM <- rbind(new_eM,c(i,ii))
							}
						}
					}
				}
				return(ftM2graphNEL(ft=cbind(V[new_eM[,1]],V[new_eM[,2]]), V=V, edgemode = "undirected"))
			}
		})



## Create a complete graph by just passing by the node names 
completegraph <- function (nodeNames) {
	
	nodeNames <- as.character(nodeNames)

	p <- length(nodeNames)
	
	adjM <- matrix(rep(1,p^2),ncol = p)
	diag(adjM) <- 0
	
	rownames(adjM) <- nodeNames
	colnames(adjM) <- nodeNames
	
	G <- new("graphAM", adjMat = adjM)
}



newgraph <- function(nodeNames, mat, weights = NULL, directed = FALSE, isAdjacency = FALSE,...) {

	nodeNames <- as.character(nodeNames)
	if(!is.matrix(mat)) {
		## mat is probably only a vector with two elements
		if(length(mat) == 2){
			mat <- matrix(mat, ncol = 2)			
		}else {
			stop("[newgraph] mat should be a matrix")
		}
	}
	edgemode <- ifelse(directed,"directed","undirected") 

	if(isAdjacency) {
		## assume it's an adjacency matrix	
		rownames(mat) <- nodeNames
		colnames(mat) <- nodeNames	
		
		G <- new("graphAM", adjMat = mat, edgemode = edgemode)
		
		if(!is.null(weights)) {
			if(!all(dim(mat)[1] == c(dim(weights)[1],dim(weights)[2]))) {
				stop("[newgraph] weights needs to be of the same dimension as mat.")
			}
			edgeDataDefaults(G, "weight") <- 1
			for(i in 1:dim(mat)[1]) {
				for(ii in 1:dim(mat)[2]) {
					if(i == ii || mat[i,ii] == 0) {
						next
					}
					edgeData(G, from = nodeNames[i], to = nodeNames[ii], attr = "weight") <- weights[i,ii]
				}
			}
		}
		
	}else {

		## assume it's a from to matrix
		if(is.numeric(mat) && !is.numeric(nodeNames)){
			mat <- cbind(nodeNames[mat[,1]],nodeNames[mat[,2]])
		}
		
		## add missing edges if directed
		n <- dim(mat)[1]
		if(!directed) {
			if(n != 1) { 
				for(i in 1:n) {
					isFound <- FALSE
					for(ii in (i+1):n) {
						if(all(mat[ii,] == mat[i,c(2,1)])) {
							isFound <- TRUE
							break
						}				
					}
					if(!isFound) {
						mat <- rbind(mat,mat[i,c(2,1)])
						weights <- c(weights,weights[i])
					}
				}
			}else { ## only one element
				mat <- rbind(mat,mat[1,c(2,1)])
				weights <- rep(weights,2)
			}
		}
	
		edL <- vector("list", length = length(nodeNames))
		names(edL) <- nodeNames
		for(i in 1:length(nodeNames)) {
			edL[[i]] <- mat[mat[,1] == nodeNames[i],2]
		}

		G <- new("graphNEL", nodes = nodeNames, edgeL=edL, edgemode = edgemode)			
		edgeDataDefaults(G, "weight") <- 1
		if(!is.null(weights)) {
			if(length(weights) != dim(mat)[1]) {
				stop("[newgraph] length of weights argument does not match with the dimension of the from-to-matrix.")
			}
			weightL <- vector("list", length = length(nodeNames))
			names(weightL) <- nodeNames
			for(i in 1:dim(mat)[1]) {
				edgeData(G, from = mat[i,1], to = mat[i,2], attr = "weight") <- weights[i]
			}
		}
		
	}	
	return(G)
}
