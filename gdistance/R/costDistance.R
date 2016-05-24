# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

#TODO check if coordinate systems are equal.
#TODO check if bounding box of coordinates falls inside bb of transition

setGeneric("costDistance", function(x, fromCoords, toCoords) standardGeneric("costDistance"))

setMethod("costDistance", signature(x = "TransitionLayer", fromCoords = "Coords", toCoords = "Coords"), def = function(x, fromCoords, toCoords)
	{
		return(.cd(x, fromCoords, toCoords))
	}
)

.cd <- function(x, fromCoords, toCoords)
{
	fromCoords <- .coordsToMatrix(fromCoords)
	toCoords <- .coordsToMatrix(toCoords)

	fromCells <- cellFromXY(x, fromCoords)
	toCells <- cellFromXY(x, toCoords)
	
	if(!all(!is.na(fromCells))){
		warning("some coordinates not found and omitted")
		fromCells <- fromCells[!is.na(fromCells)]
	}
	
	if(!all(!is.na(toCells))){
		warning("some coordinates not found and omitted")
		toCells <- toCells[!is.na(toCells)]
	}
	
	costDist <- matrix(NA, nrow=length(fromCoords[,1]),ncol=length(toCoords[,1]))
	rownames(costDist) <- rownames(fromCoords)
	colnames(costDist) <- rownames(toCoords)
	y <- transitionMatrix(x)
	if(isSymmetric(y)) {m <- "undirected"} else{m <- "directed"}
	adjacencyGraph <- graph.adjacency(y, mode=m, weighted=TRUE)
	
	E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight

	uniqueFromCells <- unique(fromCells)
	uniqueToCells <- unique(toCells)		
	shortestPaths <- shortest.paths(adjacencyGraph, v=uniqueFromCells, to=uniqueToCells, mode="out", algorithm="dijkstra")

	index1 <- match(fromCells,uniqueFromCells)
	index2 <- match(toCells,uniqueToCells)
	costDist[] <- shortestPaths[index1,index2]
	return(costDist)
}

setMethod("costDistance", signature(x = "TransitionLayer", fromCoords = "Coords", toCoords = "missing"), def = function(x, fromCoords)
	{
		return(.cd2(x, fromCoords))
	}
)

.cd2 <- function(x, fromCoords)
{
		fromCoords <- .coordsToMatrix(fromCoords)
		fromCells <- cellFromXY(x, fromCoords)

		if(!all(!is.na(fromCells))){
			warning("some coordinates not found and omitted")
			fromCells <- fromCells[!is.na(fromCells)]
		}
		
		costDist <- matrix(NA, nrow=length(fromCoords[,1]),ncol=length(fromCoords[,1]))
		rownames(costDist) <- rownames(fromCoords)
		colnames(costDist) <- rownames(fromCoords)
		
		if(isSymmetric(transitionMatrix(x))) {m <- "undirected"} else{m <- "directed"}
		adjacencyGraph <- graph.adjacency(transitionMatrix(x), mode=m, weighted=TRUE)
		
		E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight

		uniqueFromCells <- unique(fromCells)
		shortestPaths <- shortest.paths(adjacencyGraph, v=uniqueFromCells, to=uniqueFromCells, mode="out")
		index <- match(fromCells,uniqueFromCells)
		costDist[] <- shortestPaths[index,index]
		if(m=="undirected") {costDist <- as.dist(costDist)}
		return(costDist)
}