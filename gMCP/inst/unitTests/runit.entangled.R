
test.graphTest <- function() {
	G1 <- rbind(c(0,0.5,0.5,0,0), c(0,0,1,0,0),
			c(0, 0, 0, 1-0.01, 0.01), c(0, 1, 0, 0, 0),
			c(0, 0, 0, 0, 0))
	G2 <- rbind(c(0,0,1,0,0), c(0.5,0,0.5,0,0),
			c(0, 0, 0, 0.01, 1-0.01), c(0, 0, 0, 0, 0),
			c(1, 0, 0, 0, 0))
	weights <- rbind(c(1, 0, 0, 0, 0), c(0, 1, 0, 0, 0))
	pvals <- c(0.024, 0.05, 0.01, 0.003, 0.009)
	out <- graphTest(pvals, weights, alpha=c(0.0125, 0.0125), G=list(G1, G2), verbose = TRUE)
  graph1 <- matrix2graph(G1, weights[1,])
  graph2 <- matrix2graph(G2, weights[2,])
	graph <- new("entangledMCP", subgraphs=list(graph1, graph2), weights=c(0.5,0.5))
	out_from_objects <- graphTest(pvals, alpha=0.025, graph=graph, verbose = TRUE)
	checkTrue(all.equal(out, out_from_objects))
}

test.bonferroni <- function() {
  graph <- Entangled1Maurer2012()
  pvalues <- c(0.5, 0.01, 0.02, 0.01, 0.01)
  result <- gMCP(graph, pvalues, test="Bonferroni", alpha=0.05)
  checkTrue(all(result@rejected==c(FALSE, TRUE, TRUE, FALSE, TRUE)))
  pvalues <- c(0.04, 0.01, 0.02, 0.01, 0.01)
  result <- gMCP(graph, pvalues, test="Bonferroni", alpha=0.05)
  checkTrue(all(result@rejected))
}

test.parametric <- function() {
  graph <- Entangled1Maurer2012()
  pvalues <- c(0.5, 0.01, 0.02, 0.01, 0.01)
  result <- gMCP(graph, pvalues, correlation=diag(5), test="parametric" , alpha=0.05)
  
  
}

test.simes <- function() {
  
}

test.power <- function() {
  
}