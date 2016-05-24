#' Joins two graphMCP objects
#' 
#' Creates a new graphMCP object by joining two given graphMCP objects.
#' 
#' If \code{graph1} and \code{graph2} have duplicates in the node names, the
#' nodes of the second graph will be renamed.
#' 
#' If and only if the sum of the weights of graph1 and graph2 exceeds 1, the
#' weights are scaled so that the sum equals 1.
#' 
#' A description attribute of either graph will be discarded.
#' 
#' @param graph1 A graph of class \code{\link{graphMCP}}.
#' @param graph2 A graph of class \code{\link{graphMCP}}.
#' @param xOffset A numeric specifying an offset (on the x-axis) for placing
#' the nodes and edge labels of the second graph.
#' @param yOffset A numeric specifying an offset (on the y-axis) for placing
#' the nodes and edge labels of the second graph.
#' @return A graphMCP object that represents a graph that consists of the two
#' given graphs.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{graphMCP}}
#' @keywords graphs
#' @examples
#' 
#' 
#' g1 <- BonferroniHolm(2)
#' g2 <- BonferroniHolm(3)
#' 
#' joinGraphs(g1, g2)
#' 
#' 
#' @export joinGraphs
#' 
joinGraphs <- function(graph1, graph2, xOffset=0, yOffset=200) {
	m1 <- graph2matrix(graph1)
	m2 <- graph2matrix(graph2)
	m <- bdiagNA(m1,m2)
	m[is.na(m)] <- 0
	nNames <- c(getNodes(graph1), getNodes(graph2))
	d <- duplicated(nNames)
	if(any(d)) {
		warning(paste(c("The two graphs have the following identical nodes: ", paste(nNames[d], collapse=", "), ". The nodes of the second graph will be renamed."), sep=""))
		nodes2 <- getNodes(graph2)
		i <- 1
		for (x in nNames[d]) {
			while (any(nNames==paste("H",i, sep=""))) {
				i <- i + 1
			}
			nodes2[nodes2==x] <- paste("H",i, sep="")
			i <- i + 1
		}
		nNames <- c(getNodes(graph1), nodes2)
	}
	rownames(m) <- nNames
	colnames(m) <- nNames
	graph <- matrix2graph(m)	
	weights <- c(getWeights(graph1), getWeights(graph2))
	if (sum(weights)>1) {
		weights <- weights / sum(weights)
	}
	graph <- setWeights(graph, weights=weights)
	nodeX <- c(getXCoordinates(graph1), getXCoordinates(graph2) + xOffset) 
	nodeY <- c(getYCoordinates(graph1), getYCoordinates(graph2) + yOffset) 
	names(nodeX) <- nNames
	names(nodeY) <- nNames
	graph@nodeAttr$X <- nodeX
	graph@nodeAttr$Y <- nodeY
	return(graph)
}

#' Get a subgraph
#' 
#' Given a set of nodes and a graph this function creates the subgraph
#' containing only the specified nodes.
#' 
#' 
#' @param graph A graph of class \code{\link{graphMCP}}.
#' @param subset A logical or character vector specifying the nodes in the
#' subgraph.
#' @return A subgraph containing only the specified nodes.
#' @author Kornelius Rohmeyer \email{rohmeyer@@small-projects.de}
#' @seealso \code{\link{graphMCP}}
#' @keywords print graphs
#' @examples
#' 
#' 
#' graph <- improvedParallelGatekeeping()
#' subgraph(graph, c(TRUE, FALSE, TRUE, FALSE))
#' subgraph(graph, c("H1", "H3"))
#' 
#' 
#' @export subgraph
#' 
subgraph <- function(graph, subset) {
	if (is.logical(subset)) {
		subset <- getNodes(graph)[subset]	
	}
	m <- graph@m
	w <- getWeights(graph)
	rejected <- getRejected(graph)
	if (is.character(subset)) {
		subGraph <- matrix2graph(m[subset,subset], w[subset])
		setRejected(subGraph, getNodes(subGraph)) <- rejected[subset]
		return(subGraph)
	} else {
		stop("The parameter subset must be either a logical or character vector.")
	}
}