test.graphMCP <- function() {
	
	hypotheses <- LETTERS[1:4]
	m <- matrix(runif(16), nrow=4)	
	gR <- new("graphMCP", m=m, weights=c(0.1, 0.1, 0.1, 0))
	
	
	bhG5 <- BonferroniHolm(5)
	
	checkEquals(getWeights(bhG5),
				structure(c(0.2, 0.2, 0.2, 0.2, 0.2),
						.Names = c("H1", "H2", "H3", "H4", "H5")))
	
}

test.bonferroniHolm <- function() {	
	bhG3 <- BonferroniHolm(3)				   
	result <- gMCP(bhG3, pvalues=c(0.1, 0.3, 0.7), alpha=0.6)
	checkEquals(gMCP:::getRejected(result@graphs[[3]]),
			structure(c(TRUE, TRUE, FALSE), .Names = c("H1", "H2", "H3")))
}

test.gMCP <- function() {
	bhG3 <- BonferroniHolm(3)
	pvalues <- c(0.1, 0.2, 0.3)
	names(pvalues) <- getNodes(bhG3)
	checkTrue(gMCP:::canBeRejected(bhG3, "H1", alpha=0.6, pvalues)) 
	checkTrue(gMCP:::canBeRejected(bhG3, "H2", alpha=0.6, pvalues)) 
	checkTrue(!gMCP:::canBeRejected(bhG3, "H3", alpha=0.6, pvalues)) 
	checkException(gMCP(bhG3, 0, alpha=0.6))
	checkException(gMCP(bhG3, rep(0,6), alpha=0.6))
	graph <- matrix2graph(matrix(c(0,0,0,0), nrow=2))
	checkEquals(unname(getRejected(gMCP(graph, pvalues=c(0.0, 0.0)))),
			c(TRUE, TRUE))
}

test.adjPValues <- function() {
	adjPValues <- gMCP:::adjPValues(BonferroniHolm(3), c(0.02,0.055,0.012))@adjPValues
	checkEquals(adjPValues, 
			structure(c(0.04, 0.055, 0.036), .Names = c("H1", "H2", "H3")))
}

test.gMCPBretzEtAl <- function() {
	graph <- BretzEtAl2011()
	pvalues <- c(0.1, 0.008, 0.005, 0.15, 0.04, 0.006)
	result <- gMCP(graph, pvalues)
	last <- result@graphs[[4]]	
	checkEquals(unname(unlist(nodeAttr(last, getNodes(last), "rejected"))),
			c(FALSE, TRUE, TRUE, FALSE, FALSE, TRUE))
}

test.gMCPParallelGatekeeping <- function() {
	graph <- parallelGatekeeping()
	graph <- rejectNode(graph, "H1")
	
	graph <- rejectNode(graph, "H4")

}

test.only.no.error <- function() {
	graph <- BonferroniHolm(3)
	pvalues <- c(0.1, 0.2, 0.3)
	graph2latex(graph)
	gMCPReport(graph)
	gMCPReport(gMCP(graph, pvalues))
}

test.graph2matrix <- function() {
	m <- matrix(runif(9)/3, nrow=3)
	diag(m) <- c(0, 0, 0)
	rownames(m) <- LETTERS[1:3]
	colnames(m) <- LETTERS[1:3]
	m2 <- graph2matrix(matrix2graph(m))
	checkEquals(m, m2)
}