randomMCPGraph <- function(V=letters[1:10], M=1:4, p=0.2) {	
	#TODO create random graph
	g <- randomGraph(V, M, p)
	edgeAttrDefaults(g, "epsilon") <- list(0)
	class(g) <- "graphMCP"
	defaultProps <- list(nodeWeight=0, rejected=FALSE)
	nodeAttrData <- new("attrData", defaults=defaultProps)
	alpha <- 0.05
	attrDataItem(nodeAttrData, x=V, attr="nodeWeight") <- alpha/length(V)	
	g@nodeAttr <- nodeAttrData	
	g@graphData[[1]] <- "directed"
	#edgeAttrDefaults(g, "labelX") <- -100
	#edgeAttrDefaults(g, "labelY") <- -100
	for (n in getNodes(g)) {
		edgeL <- edges(g)[n][[1]]
		w <- runif(length(edgeL))
		w <- w/sum(w)
		edgeAttr(g, rep(n, length(edgeL)), edgeL, "weight") <- w
	}
	return(g)
}

isValidGraph <- function(g, alpha=0.05) {
	if (!all(TRUE==all.equal(sum(getWeights(g)), alpha))) return(paste("Sum of alpha differs from ",alpha,".",sep=""))
	for (n in getNodes(g)) {
		w <- edgeWeights(g,"d")[[1]]
		if (!all(TRUE==all.equal(sum(w),0)||TRUE==all.equal(sum(w),1))) return(paste("Sum of edges from node ",n," is ",sum(w),".",sep=""))
	}
	return(TRUE)
}

randomGraph <- function() {
	set.seed(1234)
	for (i in 1:4) {
		g <- randomMCPGraph(letters[1:10])
		checkTrue(isValidGraph(g))
		p <- runif(10)/20
		result <- gMCP(g, p)
		checkTrue(isValidGraph(result@graphs[[length(result@graphs)]]))
	}
}