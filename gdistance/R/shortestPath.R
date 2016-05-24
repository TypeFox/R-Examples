# Author: Jacob van Etten jacobvanetten@yahoo.com
# International Rice Research Institute
# Date :  January 2009, code added January 2010
# Version 1.0
# Licence GPL v3

#check if Transition and RasterLayers coincide, etc.

setGeneric("shortestPath", function(x, origin, goal, ...) standardGeneric("shortestPath"))

setMethod("shortestPath", signature(x = "TransitionLayer", origin = "Coords", goal = "Coords"), 
	def = function(x, origin, goal, output="TransitionLayer")
	{
		origin <- .coordsToMatrix(origin)
		goal <- .coordsToMatrix(goal)
		return(.shortestPath(x, origin, goal, output))		
	}
)

.shortestPath <- function(x, origin, goal, output)
{
	originCells <- cellFromXY(x, origin)
	goalCells <- cellFromXY(x, goal)
	indexOrigin <- originCells 
	indexGoal <- goalCells 
  y <- transitionMatrix(x)
	if(isSymmetric(y)) {mode <- "undirected"} else{mode <- "directed"}
	adjacencyGraph <- graph.adjacency(y, mode=mode, weighted=TRUE)
	E(adjacencyGraph)$weight <- 1/E(adjacencyGraph)$weight

	shortestPaths <- get.shortest.paths(adjacencyGraph, indexOrigin, indexGoal)$vpath
	
	if(output=="TransitionLayer")
	{
		
		result <- x
		transitionMatrix(result) <- Matrix(0, ncol=ncell(x), nrow=ncell(x))			
		for(i in 1:length(shortestPaths))
		{
			sPVector <- shortestPaths[[i]]
			adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
			adj <- rbind(adj,cbind(adj[,2], adj[,1]))
			transitionMatrix(result)[adj] <- 1/length(shortestPaths) + transitionMatrix(result)[adj]
		}
	}

	if(output=="TransitionStack")
	{
		result <- x
		transitionMatrix(result) <- Matrix(0, ncol=ncell(x), nrow=ncell(x))			
		for(i in 1:length(shortestPaths))
		{
			resultNew <- result
			sPVector <- shortestPaths[[i]] 
			adj <- cbind(sPVector[-(length(sPVector))], sPVector[-1])
			adj <- rbind(adj,cbind(adj[,2], adj[,1]))
			transitionMatrix(resultNew)[adj] <- 1/length(shortestPaths)
			result <- stack(result, resultNew)
		}
		result <- result[[2:nlayers(result)]]	
	}

	if(output=="SpatialLines")
	{
		linesList <- vector(mode="list", length=length(shortestPaths))
				
		for(i in 1:length(shortestPaths))
		{
			sPVector <- shortestPaths[[i]]
			coords <- xyFromCell(x, sPVector)
			linesList[[i]] <- Line(coords)
		}
		
  # Suggested by Sergei Petrov 
		LinesObject <- mapply(Lines, slinelist = linesList,  ID = as.character(1:length(shortestPaths)), SIMPLIFY = F)
		result <- SpatialLines(LinesObject, proj4string = CRS(projection(x)))
	}

	return(result)

}