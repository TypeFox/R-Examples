test.Simes <- function() {	
	m <- matrix(0,nr=4,nc=4)
	m[1,3] <- m[2,4] <- m[3,2] <- m[4,1] <- 1
	w <- c(1/2, 1/2, 0, 0)
	p1 <- c(0.01, 0.005, 0.01, 0.5)
	p2 <- c(0.01, 0.005, 0.015, 0.022)
	a <- 0.05
	g <- matrix2graph(m, w)
	result1 <- gMCP(g, pvalues=p1, test="Simes", alpha=a)	
	result2 <- gMCP(g, pvalues=p2, test="Simes", alpha=a)
	checkEquals(unname(result1@rejected), c(TRUE, TRUE, TRUE, FALSE))
	checkEquals(unname(result2@rejected), c(TRUE, TRUE, TRUE, TRUE))
}

checkWeights <- function(graph, pvalues) {
	# Compares the weights of the gMCP-R-code, gMCP-C-code, power-C-code and parametric-R-code
	result <- gMCP(graph,  pvalues, keepWeights=FALSE)
	rejected <- getRejected(result)
	weights <- getWeights(result)
	
	result2 <- gMCP(graph,  pvalues, useC=TRUE, keepWeights=FALSE)
	rejected2 <- getRejected(result2)
	weights2 <- getWeights(result2)
	
	checkEquals(rejected, rejected2)
	checkEquals(weights, weights2)
	
	result <- gMCP(graph,  pvalues, keepWeights=TRUE)
	rejected <- getRejected(result)
	weights <- getWeights(result)	
	
	result3 <- graphTest(pvalues=pvalues, alpha=0.05, graph=substituteEps(graph))
	m3 <- attr(result, "last.G")
	weights3 <- attr(result3, "last.alphas") / 0.05
	rejected3 <- result3!=0
	
	checkEquals(unname(rejected), unname(rejected3)) # TODO fix naming
	#checkEquals(unname(weights), weights3) TODO check why NaNs occur
}

test.checkWeights <- function() {
	graphs <- list(BonferroniHolm(5),
			parallelGatekeeping(),
			improvedParallelGatekeeping(),
			BretzEtAl2011(),
			#HungEtWang2010(),
			#HuqueAloshEtBhore2011(),
			HommelEtAl2007(),
			HommelEtAl2007Simple(),
			MaurerEtAl1995(),
			improvedFallbackI(weights=rep(1/3, 3)),
			improvedFallbackII(weights=rep(1/3, 3)),
			cycleGraph(nodes=paste("H",1:4,sep=""), weights=rep(1/4, 4)),
			fixedSequence(5),
			fallback(weights=rep(1/4, 4)),
			#generalSuccessive(weights = c(1/2, 1/2)),
			simpleSuccessiveI(),
			simpleSuccessiveII(),
			#truncatedHolm(),
			BauerEtAl2001(),
			BretzEtAl2009a(),
			BretzEtAl2009b(),
			BretzEtAl2009c()#,
			#FerberTimeDose2011(times=5, doses=3, w=1/2),
			#Ferber2011(),
			#Entangled1Maurer2012(),
			#Entangled2Maurer2012(),
	)
	for (graph in graphs) {		
		p <- gMCP:::permutations(length(getNodes(graph)))
		for (i in 1:(dim(p)[1])) {
			pvalues <- p[i,]
			pvalues[pvalues==0] <- 0.00001
			checkWeights(graph, pvalues)
		}
	}
}

test.upscale <- function() {
  g <- BonferroniHolm(5)
  r1 <- gMCP(g, pvalues=c(0.01, 0.02, 0.04, 0.04, 0.7))
  # Simple Bonferroni with empty graph:
  g2 <- matrix2graph(matrix(0, nrow=5, ncol=5))
  r2 <- gMCP(g2, pvalues=c(0.01, 0.02, 0.04, 0.04, 0.7))
  # With 'upscale=TRUE' equal to BonferroniHolm:
  r3 <- gMCP(g2, pvalues=c(0.01, 0.02, 0.04, 0.04, 0.7), upscale=TRUE)
  checkEquals(r1@rejected, r3@rejected)
  checkTrue(all(r1@rejected>=r2@rejected)) # FALSE<TRUE
}