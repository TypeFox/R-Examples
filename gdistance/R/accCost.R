# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009
# Version 1.0
# Licence GPL v3

setGeneric("accCost", function(x, fromCoords) standardGeneric("accCost"))

setMethod("accCost", signature(x = "TransitionLayer", fromCoords = "Coords"), def = function(x, fromCoords)
	{
		fromCoords <- .coordsToMatrix(fromCoords) 
		fromCells <- cellFromXY(x, fromCoords)
		if(!all(!is.na(fromCells))){
			warning("some coordinates not found and omitted")
			fromCells <- fromCells[!is.na(fromCells)]
		}
		tr <- transitionMatrix(x)
		tr <- rBind(tr,rep(0,nrow(tr)))
		tr <- cBind(tr,rep(0,nrow(tr)))
	
		startNode <- nrow(tr) #extra node to serve as origin
		adjP <- cbind(rep(startNode, times=length(fromCells)), fromCells)
	
		tr[adjP] <- Inf
	
		adjacencyGraph <- graph.adjacency(tr, mode="directed", weighted=TRUE)
		E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight		
	
		shortestPaths <- shortest.paths(adjacencyGraph, v=startNode)[-startNode]

		result <- as(x, "RasterLayer")
		result <- setValues(result, shortestPaths)	
		return(result)
	}
)

